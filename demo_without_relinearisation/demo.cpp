#include "base.hpp"
#include <vector>
#include <assert.h>
#include <NTL/ZZ.h>

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
};

struct Env 
{
    NTL::ZZ q;
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
            NTL::SetCoeff(out, i, NTL::coeff(out, i) / q +1);
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
    // sample sx with coefficients in {1, p-1} uniformly at random
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

void encrypt(Cipher &out,
             NTL::ZZX const msg,
             PublicKey const &pubKey,
             Env &env,
             PolyRing const& ring)
{
    //! try to implement your encryption function.
    // according to the fv.pdf given by hong
    // ct = (b*u+e+[q/t]*m, a*u+e)
    NTL::ZZX e, u, temp;
    ring.sampleGaussian(env.stdev, e);
    ring.sampleSmall(u);
    temp = ring.mulMod(pubKey.bx, u);
    temp = ring.addMod(temp, e);
    temp = ring.addMod(temp, (env.q/env.t)*msg);
    out.parts.push_back(temp);
    ring.sampleGaussian(env.stdev, e);
    temp = ring.mulMod(pubKey.ax, u);
    temp = ring.addMod(temp, e);
    out.parts.push_back(temp);
}

void decrypt(NTL::ZZX &out,
             Cipher const& ct,
             PrivateKey const &priKey,
             Env &env,
             PolyRing const& ring)
{
    //! try to implement your decryption function.
    
    NTL::ZZX psx = NTL::ZZX(1);
    assert(ct.parts.size() >= 2);
    // out = ring.addMod(ct.parts[0], ring.mulMod(ct.parts[1], priKey.sx));
    out = ct.parts[0];
    for (int i=1; i<ct.parts.size(); ++i)
    {
        psx = ring.mulMod(psx, priKey.sx);
        out = ring.addMod(out, ring.mulMod(psx, ct.parts[i]));
    }

    if(ct.parts.size() == 2)
    {
        ring.normalize(out, true);
        out = env.t * out;
        divRounding(out, env.q);
        modZZXByZZ(out, env.t);
    }
    else
    {
        // std::cout << "The cipher text:\n";
        // std::cout << out << std::endl;
        ring.normalize(out, true);
        out = (env.t) * (env.t) * out;
        divRounding(out, env.q * env.q);
        modZZXByZZ(out, env.t);
    }
    
}

void addCtxt(Cipher &out,
             Cipher const& ct0,
             Cipher const& ct1,
             PublicKey const& pubKey,
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
             PublicKey const& pubKey,
            //  RelinKey const& reliKey,
             PolyRing const& ring)
{
    //! try to implement your homomorphic multiplication function.
    //! If you can implemet the relinearization, it would be a plus.
    out.parts = {};
    for (int i=0; i<ct0.parts.size(); ++i)
    {
        for(int j=0; j<ct1.parts.size(); ++j)
        {
            if(i+j >= out.parts.size())
            {
                out.parts.push_back(ring.mulMod(ct0.parts[i], ct1.parts[j]));
            }
            else
            {
                out.parts[i+j] = ring.addMod(out.parts[i+j], ring.mulMod(ct0.parts[i], ct1.parts[j]));
            }
        }
    }
}

int main() {
    NTL::ZZ Q = NTL::power_ZZ(2, 120);
    NTL::ZZ T = NTL::power_ZZ(2, 90);
    const long N = 1024;
    const long hwt = 64;                // the real degree of private key
    const double stdev = 3.2;           // [-B, B], in this case, B is 8*stdev
    PolyRing ring(Q, N);                // Z_q[X]/(X^N + 1)
    Env env ={Q, T, N, hwt, stdev};

    PrivateKey prikey;                  // the coefficents of private key are {-1, 1}, padding with 0
    PublicKey pubKey;
    genPrivateKey(prikey, hwt, ring);
    genPublicKey(pubKey, stdev, prikey, ring);

    //! Note that the msg bound is [0, 2^60], considering the approximation error
    // I suggest that the msg bound is [0, 2^50]
    NTL::ZZX msg0, msg1, ptext;
    Cipher ct0, ct1, cta, ctm;

    // the initiation of msg0 and msg1
    for(int i=0; i<10; ++i){
        NTL::SetCoeff(msg0, i, 2*i);
        NTL::SetCoeff(msg1, i, i);
    }
    NTL::SetCoeff(msg0, 0, NTL::power_ZZ(2,20));
    NTL::SetCoeff(msg1, 0, NTL::power_ZZ(2,20));
    
    std::cout << ring.addMod(msg0, msg1) << std::endl;
    std::cout << ring.mulMod(msg0, msg1) << std::endl;

    encrypt(ct0, msg0, pubKey, env, ring);
    encrypt(ct1, msg1, pubKey, env, ring);

    addCtxt(cta, ct0, ct1, pubKey, ring);
    decrypt(ptext, cta, prikey, env, ring);


    mulCtxt(ctm, ct0, ct1, pubKey, ring);
    decrypt(ptext, ctm, prikey, env, ring);
    
    return 0;
}
