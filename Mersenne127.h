#pragma once

#include <x86intrin.h>
#include <gmp.h>
#include <iostream>

class Mersenne127
{
	mpz_t m_value;
public:
	Mersenne127(const char * value = NULL);
	Mersenne127(const mpz_t value);
	Mersenne127(const Mersenne127 & other);
	Mersenne127(int);
	~Mersenne127();

    const Mersenne127 & operator = (const Mersenne127 & other);
    bool operator == (const Mersenne127 & other) const;
    bool operator != (const Mersenne127 & other) const;

    Mersenne127 operator + (const Mersenne127 & rha) const;
    Mersenne127 operator - (const Mersenne127 & rha) const;
    Mersenne127 operator * (const Mersenne127 & rha) const;
    Mersenne127 operator / (const Mersenne127 & rha) const;

    const Mersenne127 & operator += (const Mersenne127 & rha);
    const Mersenne127 & operator *= (const Mersenne127 & rha);

    const mpz_t * get_mpz_t() const { return &m_value; }
    void set_mpz_t(const mpz_t * value) { mpz_set(m_value, *value); }

    friend std::ostream& operator << (std::ostream& s, const Mersenne127 & m127);
    friend std::istream& operator >> (std::istream& s, Mersenne127 & m127);

    static const char M127[];
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

class ZpMersenne127Element {

//private:
public: //TODO return to private after tesing

    //we use gcc 5.4 supported uint 128 bit type
    __uint128_t elem;
    static const __uint128_t p;
    //the prime is 2^127-1

public:

    ZpMersenne127Element(){elem = 0;};
    ZpMersenne127Element(__uint128_t e)
    {
        this->elem = e;
        if(this->elem>=p){

            this->elem-= p;
        }
    }

    inline ZpMersenne127Element& operator=(const ZpMersenne127Element& other){
        elem = other.elem;
        return *this;
    }
    inline bool operator!=(const ZpMersenne127Element& other)

    { return !(other.elem == elem); };

    inline bool operator==(const ZpMersenne127Element& other)

    { return other.elem == elem; };

    ZpMersenne127Element operator+(const ZpMersenne127Element& f2)
    {
        ZpMersenne127Element answer;
        answer.elem = (elem + f2.elem);
        if(answer.elem>=p)
            answer.elem-=p;
        return answer;
    }

    ZpMersenne127Element operator-(const ZpMersenne127Element& f2)
    {
        ZpMersenne127Element answer;
        __int128_t temp =  (__int128_t)elem - (__int128_t)f2.elem;
        if(temp<0){
            answer.elem = temp + p;
        }
        else{
            answer.elem = temp;
        }
        return answer;
    }

    ZpMersenne127Element operator/(const ZpMersenne127Element& f2)
    {
        mp_limb_t* me = (mp_limb_t *) &elem;

        auto f2Eleme = f2.elem;
        mp_limb_t* elemOff2 = (mp_limb_t *) &(f2Eleme); //keep same naming convention as other ZpMersenne classes
        mp_limb_t *d = (mp_limb_t *) &p; //d is actually p
        mp_limb_t result_1[2]; //result is used a few times. we do not allow override in low-level, so more vars.
        mp_limb_t result_2[4]; //result of mult is 256 bits (limb is 64 bit)
        __uint128_t res;

        mp_limb_t tp[16]; //scratch space for invert
        mpn_sec_invert (result_1, elemOff2, d, 2, 256, tp);
        mpn_mul (result_2, result_1, 2,me, 2);

        mp_limb_t q[4];
        mpn_tdiv_qr (q, 				//quotent - not used
                     (mp_limb_t *)&res, //remainder
                     0, 				//nuat be 0
                     result_2,			//
                     4,
                     d,					//mod divisor (d == p)
                     2);

        return ZpMersenne127Element(res);
    }

    ZpMersenne127Element operator*(const ZpMersenne127Element& f2)
    {
        unsigned long long m64[4];

        //do four mults
        unsigned long long* me = (unsigned long long*)&elem;
        unsigned long long* other = (unsigned long long*)&f2.elem;

        unsigned long long high00;
        unsigned long low00 = _mulx_u64(me[0], other[0], &high00);

        unsigned long long high01;
        unsigned long low01 = _mulx_u64(me[0], other[1], &high01);

        unsigned long long high10;
        unsigned long low10 = _mulx_u64(me[1], other[0], &high10);

        unsigned long long high11;
        unsigned long low11 = _mulx_u64(me[1], other[1], &high11);

        m64[0] = low00;

        unsigned char c1 = 0, c2 = 0;
        c1 = _addcarry_u64(c1, high00, low01, &m64[1]);
        c2 = _addcarry_u64(c2, m64[1], low10, &m64[1]);
        //m64[1] = high00+low01+low10;
        //m64[2] = high01 + high10 + low11 + c1 +c2;
        //c1=c2=0;
        c1 = _addcarry_u64(c1, high01, high10, &m64[2]);
        c2 = _addcarry_u64(c2, m64[2], low11, &m64[2]);

        m64[3] = high11+c1+c2;

        __uint128_t *m128 = (__uint128_t *) &m64;

        // mpn_mul ( (mp_limb_t *)m128, (mp_limb_t *)&elem, 2,	(mp_limb_t *)&(f2.elem), 2);
        __uint128_t res, low, low128bit , hign,highShift1;

        low = (m128[0] & p);
        low128bit = (m128[0]>>127);
        highShift1 = (m128[1]<<1);
        res = low + low128bit + highShift1;

        return ZpMersenne127Element(res);
    }

    ZpMersenne127Element& operator+=(const ZpMersenne127Element& f2)
    {
        elem = (elem + f2.elem);

        if(elem>=p)
            elem-=p;

        return *this;
    }

    ZpMersenne127Element& operator*=(const ZpMersenne127Element& f2)
    {
        return *this = *this * f2;
    }

    void get_mpz_t(mpz_t value) const
    {
        mpz_import(value, 1, -1, 16, 0, 0, (void*)&elem);
        return;
    }

    void set_mpz_t(const mpz_t value)
    {
        mpz_export((void*)&elem, NULL, -1, 16, 0, 0, value);
        return;
    }

    static void get_mpz_t_p(mpz_t value)
    {
        mpz_import(value, 1, -1, 16, 0, 0, (void*)&p);
        return;
    }
};

inline std::ostream& operator<<(std::ostream& s, const ZpMersenne127Element & a){
	const u_int64_t * abytes = (const u_int64_t *)&a;
	return s << abytes[1] << " " << abytes[0] << std::endl;
};

