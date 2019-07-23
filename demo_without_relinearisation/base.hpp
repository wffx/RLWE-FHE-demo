#pragma once
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <cmath>
#include <numeric>
#include <vector>
#include <assert.h>

#define LONG_MAX 1000000000
//! Class for arithemtic over the ring Z[X] / (X^N + 1) for N is 2-power number.
class PhiMPolyBase {
public:
    PhiMPolyBase(long N) : N_(N) {
        assert(N > 0);
        assert(is_two_pow(N) && "The polynomial degree must be 2-power number");
        phiMX_.SetLength(N_ + 1);
        //phiMX_ = X^N + X^0
        NTL::SetCoeff(phiMX_, 0, 1);
        NTL::SetCoeff(phiMX_, N_, 1);
    }

    ~PhiMPolyBase() {}

    PhiMPolyBase(PhiMPolyBase const& oth) 
        : N_(oth.N_),
          phiMX_(oth.phiMX_) { }

    PhiMPolyBase& operator=(PhiMPolyBase const& oth) {
        PhiMPolyBase copy(oth);
        swap(*this, copy);
        return *this;
    }

    void swap(PhiMPolyBase &a, PhiMPolyBase &b) const {
        std::swap(a.N_, b.N_);
        std::swap(a.phiMX_, b.phiMX_);
    }

    NTL::ZZX addMod(NTL::ZZX const& a, NTL::ZZX const& b) const {
        assert(NTL::deg(a) <= N_);
        assert(NTL::deg(b) <= N_);
        NTL::ZZX ans;
        NTL::add(ans, a, b);
        return ans;
    }

    NTL::ZZX subMod(NTL::ZZX const& a, NTL::ZZX const& b) const {
        assert(NTL::deg(a) <= N_);
        assert(NTL::deg(b) <= N_);
        NTL::ZZX ans;
        NTL::sub(ans, a, b);
        return ans;
    }

    NTL::ZZX mulMod(NTL::ZZX const& a, NTL::ZZX const& b) const {
        assert(NTL::deg(a) <= N_);
        assert(NTL::deg(b) <= N_);
        NTL::ZZX ans;
        NTL::MulMod(ans, a, b, phiMX_);
        return ans;
    }

    NTL::ZZX mulMod(NTL::ZZX const& a, NTL::ZZ const& b) const {
        assert(NTL::deg(a) <= N_);
        NTL::ZZX ans;
        ans.SetLength(a.rep.length());
        for (long i = 0; i < a.rep.length(); ++i)
            NTL::SetCoeff(ans, i, NTL::coeff(a, i) * b);
        return ans;
    }

    long getN() const { return N_; }

private:
    bool is_two_pow(long v) const {
        return (v & -v) == v;
    }

private:
    long N_;
    NTL::ZZX phiMX_;
};

//! Class for polynomial ring Z_p[X]/(X^N + 1) for positive value p and N.
class PolyRing  {
public:
    explicit PolyRing(NTL::ZZ const& p, long N) 
        : p_(p),
          phalf_(p >> 1),
          base_(N) {
        assert(p >= 2 && "The moduli should not be less than 2.");
    }
    
    PolyRing(long const& p, long N) : PolyRing(NTL::to_ZZ(p), N) {}

    const NTL::ZZ getP() const { return p_; }

    const long getN() const { return base_.getN(); }

    ~PolyRing() {

    }

    PolyRing(PolyRing const& oth) : p_(oth.p_), phalf_(oth.phalf_), base_(oth.base_) {}

    void swap(PolyRing &a, PolyRing &b) {
        std::swap(a.p_, b.p_);
        std::swap(a.phalf_, b.phalf_);
        std::swap(a.base_, b.base_);
    }

    PolyRing& operator=(const PolyRing &oth) {
        PolyRing copy(oth);
        swap(*this, copy);
        return *this;
    }

    NTL::ZZX addMod(NTL::ZZX const& a, NTL::ZZX const& b) const {
        NTL::ZZX ans = base_.addMod(a, b);
        normalize(ans);
        return ans;
    }

    NTL::ZZX subMod(NTL::ZZX const& a, NTL::ZZX const& b) const {
        NTL::ZZX ans = base_.subMod(a, b);
        normalize(ans);
        return ans;
    }

    NTL::ZZX mulMod(NTL::ZZX const& a, NTL::ZZX const& b) const {
        NTL::ZZX ans = base_.mulMod(a, b);
        normalize(ans);
        return ans;
    }

    NTL::ZZX mulMod(NTL::ZZX const& a, NTL::ZZ const& b) const {
        NTL::ZZX ans = base_.mulMod(a, b);
        normalize(ans);
        return ans;
    }

    NTL::ZZX sampleSmall() const {
        NTL::ZZX ans;
        sampleSmall(ans);
        return ans;
    }

    //! sample polynomial with coefficients in {-1, 0, 1} uniformly at random.
    void sampleSmall(NTL::ZZX &ans) const {
        const long n = base_.getN();
        ans.SetLength(n);
        for (long i = 0; i < n; ++i) {
            long r = NTL::RandomBits_long(2); // {0, 1, 2, 3}
            if(r & 1)
                NTL::SetCoeff(ans, i, 0);
            else
                NTL::SetCoeff(ans, i, (r & 2) - 1); // {-1, 1}
        }
        normalize(ans);
    }

    NTL::ZZX sampleUniform() const {
        NTL::ZZX ans;
        sampleUniform(ans);
        return ans;
    }

    NTL::ZZX sampleHWT(long k) const {
        NTL::ZZX ans;
        sampleHWT(k, ans);
        return ans;
    }

    // sample polynomial with K coefficients in {p-1, 1}, others padding with 0 
    void sampleHWT(long k, NTL::ZZX &ans) const {
        const long n = base_.getN();
        assert(k > 0 && k <= n);
        std::vector<long> pool(k);
        std::iota(pool.begin(), pool.end(), 0);
        for (long i = k; i < n; ++i) {
            long j = NTL::RandomBnd(LONG_MAX) % i;
            if (j < k)
                pool.at(j) = i;
        }

        ans.rep.kill();
        ans.SetLength(n);
        for (long pos : pool) {
            long r = NTL::RandomBits_long(2) & 2; // {0, 2}
            NTL::SetCoeff(ans, pos, r - 1); // {-1, 1}
        }
        normalize(ans);
    }

    void sampleUniform(NTL::ZZX &ans) const {
        const long n = base_.getN();
        ans.SetLength(n);
        for (long i = 0; i < n; ++i)
            NTL::SetCoeff(ans, i, NTL::RandomBnd(p_));
    }

    void sampleGaussian(double stdev, NTL::ZZX &ans) const {
        assert(stdev > 0.);
        static double const Pi= 4.0 * std::atan(1.0);  // Pi=3.1415..
        static double const bignum = (double) LONG_MAX; // convert to double

        const long n = base_.getN();
        ans.rep.kill();
        ans.SetLength(n);
        for (long i = 0; i < n; i += 2) {
            // r1, r2 are "uniform in (0,1)"
            double r1 = (1 + NTL::RandomBnd(LONG_MAX)) / (bignum + 1);
            double r2 = (1 + NTL::RandomBnd(LONG_MAX)) / (bignum + 1);
            double theta = 2 * Pi * r1;
            double rr = std::sqrt(-2.0 * std::log(r2)) * stdev;
            if (rr > 8 * stdev) // sanity-check, trancate at 8 standard deviations
                rr = 8 * stdev;

            // Generate two Gaussians RV's
            double g0 = rr * std::cos(theta);
            double g1 = rr * std::sin(theta);
            NTL::SetCoeff(ans, i, (long) g0);
            if (i + 1 < n)
                NTL::SetCoeff(ans, i, (long) g1);
        }

        normalize(ans);
    }

    NTL::ZZX sampleGaussian(double stdev) const {
        NTL::ZZX ans;
        sampleGaussian(stdev, ans);
        return ans;
    }

    //! correct the coefficients of the polynomial to Z_p
    //! When negative is true, reduce to [-p/2, p/2)
    void normalize(NTL::ZZX &ply, bool negative = false) const {
        const long n = ply.rep.length();
        for (long i = 0; i < n; ++i) {
            NTL::SetCoeff(ply, i, NTL::coeff(ply, i) % p_);
            if (NTL::coeff(ply, i) < 0)
                ply[i] += p_;
        }

        if (negative) {
            for (long i = 0; i < n; ++i) {
                NTL::ZZ const& c = NTL::coeff(ply, i);
                if (c >= phalf_)
                    NTL::SetCoeff(ply, i, c - p_);
            }
        }
    }
protected:
    NTL::ZZ p_;
    NTL::ZZ phalf_;
    PhiMPolyBase base_;
};
