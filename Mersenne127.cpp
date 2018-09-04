#include "Mersenne127.h"
#include <string>
#include <memory.h>

#include "TemplateField.h"

const char Mersenne127::M127[] = "170141183460469231731687303715884105727";

class Mersenne127Helper
{
public:
	mpz_t m_M127;

	Mersenne127Helper()
	{
		mpz_init(m_M127);
		mpz_set_str(m_M127, Mersenne127::M127, 10);
	}

	~Mersenne127Helper()
	{
		mpz_clear(m_M127);
	}
};

static Mersenne127Helper the_help;

Mersenne127::Mersenne127(const char * value/* = NULL*/)
{
	mpz_init(m_value);
	if(NULL != value)
	{
		mpz_set_str(m_value, value, 10);
		mpz_mod(m_value, m_value, the_help.m_M127);
	}
}

Mersenne127::Mersenne127(const mpz_t value)
{
	mpz_init(m_value);
	mpz_mod(m_value, value, the_help.m_M127);
}

Mersenne127::Mersenne127(const Mersenne127 & other)
{
	mpz_init(m_value);
	mpz_mod(m_value, other.m_value, the_help.m_M127);
}

Mersenne127::Mersenne127(int i)
{
	mpz_init(m_value);
	mpz_set_si(m_value, i);
}

Mersenne127::~Mersenne127()
{
	mpz_clear(m_value);
}

const Mersenne127 & Mersenne127::operator = (const Mersenne127 & other)
{
	if(&other != this)
	{
		mpz_set(m_value, other.m_value);
	}
	return *this;
}

bool Mersenne127::operator == (const Mersenne127 & other) const
{
	return (0 == mpz_cmp(m_value, other.m_value))? true: false;
}

bool Mersenne127::operator != (const Mersenne127 & other) const
{
	return (0 != mpz_cmp(m_value, other.m_value))? true: false;
}

Mersenne127 Mersenne127::operator + (const Mersenne127 & rha) const
{
	Mersenne127 x;

	mpz_add(x.m_value, m_value, rha.m_value);
	mpz_mod(x.m_value, x.m_value, the_help.m_M127);
	return x;
}

Mersenne127 Mersenne127::operator - (const Mersenne127 & rha) const
{
	Mersenne127 x;
	mpz_sub(x.m_value, m_value, rha.m_value);
	mpz_mod(x.m_value, x.m_value, the_help.m_M127);
	return x;
}

Mersenne127 Mersenne127::operator * (const Mersenne127 & rha) const
{
	Mersenne127 x;
	mpz_mul(x.m_value, m_value, rha.m_value);
	mpz_mod(x.m_value, x.m_value, the_help.m_M127);
	return x;
}

Mersenne127 Mersenne127::operator / (const Mersenne127 & rha) const
{
	Mersenne127 x;
    //Mersenne127 inverse;
//	mpz_cdiv_q(x.m_value, m_value, rha.m_value);
//	mpz_mod(x.m_value, x.m_value, the_help.m_M127);


	mpz_invert(x.m_value, rha.m_value, the_help.m_M127);
    mpz_mul(x.m_value, x.m_value, m_value);
    mpz_mod(x.m_value, x.m_value, the_help.m_M127);



	return x;
}

const Mersenne127 & Mersenne127::operator += (const Mersenne127 & rha)
{
	mpz_add(m_value, m_value, rha.m_value);
    mpz_mod(m_value, m_value, the_help.m_M127);
	return *this;
}

const Mersenne127 & Mersenne127::operator *= (const Mersenne127 & rha)
{
	mpz_mul(m_value, m_value, rha.m_value);
    mpz_mod(m_value, m_value, the_help.m_M127);
	return *this;
}

std::ostream& operator << (std::ostream& s, const Mersenne127 & m127)
{
    return s<<m127.m_value;
};



std::istream& operator >> (std::istream& s, Mersenne127 & m127)
{
	std::string str_m127;
	s >> str_m127;
	mpz_set_str(m127.m_value, str_m127.c_str(), 10);
	return s;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

template <>
TemplateField<ZpMersenne127Element>::TemplateField(long fieldParam) {

    ZpMersenne127Element::init();

    this->elementSizeInBytes = 16;//round up to the next byte
    this->elementSizeInBits = 127;

    auto randomKey = prg.generateKey(128);
    prg.setKey(randomKey);

    m_ZERO = new ZpMersenne127Element(0);
    m_ONE = new ZpMersenne127Element(1);


}

template <>
ZpMersenne127Element TemplateField<ZpMersenne127Element>::GetElement(long b) {


    if(b == 1)
    {
        return *m_ONE;
    }
    if(b == 0)
    {
        return *m_ZERO;
    }
    else{
        ZpMersenne127Element element(b);
        return element;
    }
}

template <>
void TemplateField<ZpMersenne127Element>::elementToBytes(unsigned char* elemenetInBytes, ZpMersenne127Element& element){

    memcpy(elemenetInBytes, (byte*)(&element.elem), 16);
}

/* Not used in this present library
template <>
void TemplateField<ZpMersenne127Element>::elementVectorToByteVector(vector<ZpMersenne127Element> &elementVector, vector<byte> &byteVector){

    copy_byte_array_to_byte_vector((byte *)elementVector.data(), elementVector.size()*elementSizeInBytes, byteVector,0);
}
*/

template <>
ZpMersenne127Element TemplateField<ZpMersenne127Element>::bytesToElement(unsigned char* elemenetInBytes){

    return ZpMersenne127Element((__uint128_t)(*(__uint128_t *)elemenetInBytes));
}

__uint128_t ZpMersenne127Element::p = 0;
