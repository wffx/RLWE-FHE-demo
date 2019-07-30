#include "base.hpp"
#include <vector>

struct Cipher {
    std::vector<NTL::ZZX> parts;
};

struct PrivateKey {
    NTL::ZZX sx;
};

struct PublicKey {
    NTL::ZZX bx, ax;
};

struct RelinKey {
    //! try to implement your multiplication key.
    NTL::ZZX rlkx[12][2];
};

struct Env 
{
    NTL::ZZ Q;
    NTL::ZZ T;
    NTL::ZZ t;
    long N;
    long hwt;
    double stdev;
};

void divRounding(NTL::ZZX &out, NTL::ZZ q)
{
    NTL::ZZ r;
    for(int i=0; i<out.rep.length(); ++i)           // div by q
    {
        // NTL::SetCoeff(out, i, NTL::coeff(out, i) / NTL::ZZ(env.q));
        NTL::rem(r, NTL::coeff(out, i), NTL::ZZ(q));
        if (r < q/2)
        {
            NTL::SetCoeff(out, i, NTL::coeff(out, i) / q);
        }
        else
        {
            NTL::SetCoeff(out, i, NTL::coeff(out, i) / q + 1);
        }
    }
}

void modZZXByZZ(NTL::ZZX &out, NTL::ZZ t)
{
    for(int i=0; i<out.rep.length(); ++i)           // mod by t
    {
        NTL::SetCoeff(out, i, NTL::coeff(out, i) % t);
    }
}

void genPrivateKey(PrivateKey &priKey, 
                   const long hwt,
                   PolyRing const& ring)
{
    ring.sampleHWT(hwt, priKey.sx);
}



void genPublicKey(PublicKey &pubKey, 
                  const double stdev,
                  PrivateKey const& priKey,
                  PolyRing const& ring)
{
    // pk := (-ax * sk + e, ax)
    ring.sampleUniform(pubKey.ax);
    ring.sampleGaussian(stdev, pubKey.bx);
    pubKey.bx = ring.subMod(pubKey.bx, ring.mulMod(pubKey.ax, priKey.sx));
}

void getRelinKey(RelinKey &reliKey,
                 PrivateKey const &prikey,
                 PolyRing const& ring,
                 const Env &env)
{
    NTL::ZZX ss;
    ss = ring.mulMod(prikey.sx, prikey.sx);

    for (int i=0; i<12 ; ++i) {
        ring.sampleUniform(reliKey.rlkx[i][1]);
        ring.sampleGaussian(env.stdev, reliKey.rlkx[i][0]);
        reliKey.rlkx[i][0] = ring.subMod(reliKey.rlkx[i][0], ring.mulMod(reliKey.rlkx[i][1], prikey.sx));
        reliKey.rlkx[i][0] = ring.addMod(reliKey.rlkx[i][0], power(env.T,i)*ss);
    }
}
void encrypt(Cipher &out,
             NTL::ZZX const msg,
             PublicKey const &pubKey,
             PolyRing const& ring,
             const Env &env)
{
    //! try to implement your encryption function.
    NTL::ZZX ux, ex, temp;

    ring.sampleGaussian(env.stdev, ux);
    ring.sampleGaussian(env.stdev, ex);
    temp = ux;
    temp = ring.mulMod(temp, pubKey.bx);
    temp = ring.addMod(temp, ex);
    temp = ring.addMod(temp, env.Q/env.t*msg);
    out.parts.push_back(temp);
    ring.sampleGaussian(env.stdev, ex);
    ux = ring.mulMod(ux, pubKey.ax);
    ux = ring.addMod(ux, ex);
    out.parts.push_back(ux);
}

void decrypt(NTL::ZZX &out,
             Cipher const& ct,
             PrivateKey const &priKey,
             PolyRing const& ring,
             const Env env)
{
    //! try to implement your decryption function.
    NTL::ZZX psx = NTL::ZZX(1);
    assert(ct.parts.size() >= 2);

    out = ct.parts[0];
    for (int i=1; i<ct.parts.size(); ++i)
    {
        psx = ring.mulMod(psx, priKey.sx);
        out = ring.addMod(out, ring.mulMod(psx, ct.parts[i]));
    }

    // ring.normalize(out, true);
    divRounding(out, env.Q/env.t);
    modZZXByZZ(out, env.t);
    
    if(ct.parts.size() >2) {
        // std::cout << out << std::endl;
    }
}

void addCtxt(Cipher &out,
             Cipher const& ct0,
             Cipher const& ct1,
            //  PublicKey const& pubKey,
             PolyRing const& ring)
{
    //! try to implement your homomorphic addition function.
    out.parts = {};
    long ll;
    ll = ct0.parts.size() > ct1.parts.size() ? ct0.parts.size() : ct1.parts.size();
    out = ct0.parts.size() > ct1.parts.size() ? ct0 : ct1;
    for(int i=0; i<ll; ++i)
    {
        if (ct0.parts.size() > ct1.parts.size())
        {
           out.parts[i] = ring.addMod(out.parts[i], ct1.parts[i]);
        }
        else
        {
            out.parts[i] = ring.addMod(out.parts[i], ct0.parts[i]);
        }
    }
}

void mulCtxt(Cipher &out,
             Cipher const& ct0,
             Cipher const& ct1,
            //  PublicKey const& pubKey,
            //  RelinKey const& reliKey,
             PolyRing const& ring,
             const Env &env)
{
    //! try to implement your homomorphic multiplication function.
    //! If you can implemet the relinearization, it would be a plus.
    out.parts = {};
    NTL::ZZX temp = NTL::ZZX(0);
    
    for (int k=0; k < ct0.parts.size()+ct1.parts.size()-1; ++k) {
        for(int i=0; i<k; ++i) {
            if (i >= ct0.parts.size() || (k-i) >= ct1.parts.size()) {
                continue;
            }
            temp =ring.addMod(temp, ring.mulMod(ct0.parts[i], ct1.parts[k-i]));
        }
        ring.normalize(temp, true);
        divRounding(temp, env.Q/env.t);
        out.parts.push_back(temp);
    }
}

void reliCtxt(const RelinKey &reliKey, 
              Cipher &ct,
              const PolyRing &ring,
              const Env &env)
{
    std::vector<NTL::ZZX> ptx;
    NTL::ZZX temp, px;
    NTL::ZZ d;
    assert(ct.parts.size() == 3);

    temp = ct.parts[2];
    px = temp;
    for (int i=11; i>=0; --i) {
        
        for (int j=0; j<temp.rep.length(); ++j) {
            NTL::div(d, NTL::coeff(temp, j), NTL::power(env.T, i));
            NTL::SetCoeff(temp, j, NTL::coeff(temp, j)-d*NTL::power(env.T, i));
            NTL::SetCoeff(px, j, d);
        }
        ptx.push_back(px);
    }
    
    for(int i=0; i < ptx.size()/2; ++i) {
        px = ptx[i];
        ptx[i] = ptx[ptx.size()-i];
        ptx[ptx.size()-i] = px;
    }
    std::cout << "----" << ptx.size() << "----" << std::endl;
    temp = 0;
    for (int i=0; i<12; ++i) {
        temp = temp + power(env.T, i)*ptx[i];
    }
    if (temp == ct.parts[2]) {
        std::cout << "Correct!!!" << std::endl;
    }

    temp = 0; px = 0; 

    for(int i=0; i<12; ++i) {
        temp = temp + ring.mulMod(reliKey.rlkx[i][0], ptx[i]);
        px = px + ring.mulMod(reliKey.rlkx[i][1], ptx[i]); 
    }

    ct.parts[0] = ct.parts[0] + temp;
    ct.parts[1] = ct.parts[1] + px;

    ct.parts.pop_back();
}

int main() {
    NTL::ZZ Q = NTL::power_ZZ(2, 120);
    NTL::ZZ T = NTL::power_ZZ(2, 10);
    NTL::ZZ t = NTL::power_ZZ(2, 50);
    const long N = 1024;
    const long hwt = 64;
    const double stdev = 3.2;
    PolyRing ring(Q, N); // Z_Q[X]/(X^N + 1)
    Env env ={Q, T, t, N, hwt, stdev};

    PrivateKey priKey;
    PublicKey pubKey;
    RelinKey reliKey;
    genPrivateKey(priKey, hwt, ring);
    genPublicKey(pubKey, stdev, priKey, ring);
    getRelinKey(reliKey, priKey, ring, env);
    
    NTL::ZZX msg0, msg1, ptext;
    Cipher ct0, ct1, cta, ctm;
    
    for(int i=0; i<10; ++i){
        NTL::SetCoeff(msg0, i, 2*i);
        NTL::SetCoeff(msg1, i, i);
    }

    encrypt(ct0, msg0, pubKey, ring, env);
    encrypt(ct1, msg1, pubKey, ring, env);

    addCtxt(cta, ct0, ct1, ring);
    decrypt(ptext, cta, priKey, ring, env);
    std::cout << ptext << std::endl;

    mulCtxt(ctm, ct0, ct1, ring, env);
    std::cout << "xxxx" << std::endl;
    // reliCtxt(reliKey, ctm, ring, env);
    // std::cout << ctm.parts.size() << std::endl;
    // std::cout << ctm.parts[1] << std::endl;
    decrypt(ptext, ctm, priKey, ring, env);

    std::cout << ptext << std::endl;
    return 0;
}
