//
// Created by meital on 15/11/16.
//

#include "TemplateField.h"

#include "ZpMersenneLongElement.h"
#include "Mersenne127.h"


static const int order = -1;
static const int size = 8;
static const int endian = 0;
static const int nails = 0;

using namespace NTL;

template <>
TemplateField<ZZ_p>::TemplateField(long fieldParam) {

    this->fieldParam = fieldParam;
    this->elementSizeInBytes = NumBytes(fieldParam);//round up to the next byte
    this->elementSizeInBits = this->elementSizeInBytes*8;

    ZZ_p::init(ZZ(fieldParam));

    auto randomKey = prg.generateKey(128);
    prg.setKey(randomKey);

    m_ZERO = new ZZ_p(0);
    m_ONE = new ZZ_p(1);
}


template <>
TemplateField<ZpMersenneLongElement>::TemplateField(long fieldParam) {

    this->elementSizeInBytes = 8;//round up to the next byte
    this->elementSizeInBits = 61;

    auto randomKey = prg.generateKey(128);
    prg.setKey(randomKey);

    m_ZERO = new ZpMersenneLongElement(0);
    m_ONE = new ZpMersenneLongElement(1);
}



/*
 * The i-th field element. The ordering is arbitrary, *except* that
 * the 0-th field element must be the neutral w.r.t. addition, and the
 * 1-st field element must be the neutral w.r.t. multiplication.
 */
template <>
GF2E TemplateField<GF2E>::GetElement(long b) {

    if(b == 1)
    {
        return *GetOne();
    }
    if(b == 0)
    {
        return *GetZero();
    }
    GF2X element;

    for(int i=0; i < fieldParam; i++) {
        // set the coefficient of x^i to 1
        SetCoeff(element,i,(b >> i) & 1);
    }

    return to_GF2E(element);
}

template <>
ZpMersenneLongElement TemplateField<ZpMersenneLongElement>::GetElement(long b) {


    if(b == 1)
    {
        return *m_ONE;
    }
    if(b == 0)
    {
        return *m_ZERO;
    }
    else{
        ZpMersenneLongElement element(b);
        return element;
    }
}



template <>
ZZ_p TemplateField<ZZ_p>::GetElement(long b) {


    if(b == 1)
    {
        return *m_ONE;
    }
    if(b == 0)
    {
        return *m_ZERO;
    }
    else{
        ZZ_p element(b);
        return element;
    }
}



/**
 * the function create a field by:
 * generate the irreducible polynomial x^8 + x^4 + x^3 + x + 1 to work with
 * init the field with the newly generated polynomial
 */
template <>
TemplateField<GF2E>::TemplateField(long fieldParam) {

    this->fieldParam = fieldParam;
    this->elementSizeInBytes = fieldParam/8;
    this->elementSizeInBits = elementSizeInBytes*8;
    GF2X irreduciblePolynomial = BuildSparseIrred_GF2X(fieldParam);
    GF2E::init(irreduciblePolynomial);

    auto randomKey = prg.generateKey(128);
    prg.setKey(randomKey);

    m_ZERO = new GF2E(0);
    m_ONE = new GF2E(1);
}


template <>
void TemplateField<GF2E>::elementToBytes(unsigned char* elemenetInBytes, GF2E& element){

    BytesFromGF2X(elemenetInBytes,rep(element),fieldParam/8);
}



template <>
void TemplateField<ZpMersenneLongElement>::elementToBytes(unsigned char* elemenetInBytes, ZpMersenneLongElement& element){

    memcpy(elemenetInBytes, (byte*)(&element.elem), 8);
}

template <>
ZpMersenneLongElement TemplateField<ZpMersenneLongElement>::bytesToElement(unsigned char* elemenetInBytes){

    return ZpMersenneLongElement((unsigned long)(*(unsigned long *)elemenetInBytes));
}

template <>
ZpMersenneLongElement TemplateField<ZpMersenneLongElement>::Random() {
    unsigned long b = (prg.getRandom64() >> 3);
    return GetElement(b);
}

template <>
GF2E TemplateField<GF2E>::bytesToElement(unsigned char* elemenetInBytes){

    //first create a GF2X
    GF2X polynomialElement;

    //translate the bytes into a GF2X element
    GF2XFromBytes(polynomialElement, elemenetInBytes, fieldParam/8);

    //convert the GF2X to GF2E
    return to_GF2E(polynomialElement);
}

template <>
GF2E TemplateField<GF2E>::Random() {
    unsigned long b;
    if(elementSizeInBytes<=4)
        b = prg.getRandom32();
    else {
        b = prg.getRandom64()>>(64-elementSizeInBits);
    }
    return GetElement(b);
}

template <>
void TemplateField<ZZ_p>::elementToBytes(unsigned char* elemenetInBytes, ZZ_p& element){

    BytesFromZZ(elemenetInBytes,rep(element),elementSizeInBytes);
}


template <>
ZZ_p TemplateField<ZZ_p>::bytesToElement(unsigned char* elemenetInBytes){

    //first create a ZZ
    ZZ zz;

    //translate the bytes into a ZZ element
    ZZFromBytes(zz, elemenetInBytes, elementSizeInBytes);

    //convert the ZZ to ZZ_p
    return to_ZZ_p(zz);
}


template <>
TemplateField<Mersenne127>::TemplateField(long fieldParam)
{
    this->elementSizeInBytes = 16;//round up to the next byte
    this->elementSizeInBits = 127;

    auto randomKey = prg.generateKey(128);
    prg.setKey(randomKey);

    m_ZERO = new Mersenne127(0);
    m_ONE = new Mersenne127(1);
}

template <>
void TemplateField<Mersenne127>::elementToBytes(unsigned char* elemenetInBytes, Mersenne127& element)
{
    memset(elemenetInBytes, 0, 16);
    mpz_export((void*)elemenetInBytes, NULL, order, size, endian, nails, *element.get_mpz_t());
}

template <>
Mersenne127 TemplateField<Mersenne127>::bytesToElement(unsigned char* elemenetInBytes)
{
    mpz_t value;
    mpz_init(value);
    mpz_import(value, 2, order, size, endian, nails, (void*)elemenetInBytes);
    Mersenne127 element(value);
    mpz_clear(value);
    return element;
}

template <>
Mersenne127 TemplateField<Mersenne127>::GetElement(long b)
{
    if(b == 1)		return *m_ONE;
    if(b == 0)		return *m_ZERO;
    else
    {
        Mersenne127 element(b);
        return element;
    }
}

template <>
ZpMersenne127Element TemplateField<ZpMersenne127Element>::Random() {
	__uint128_t v = ((__uint128_t)prg.getRandom64())<<64 | (__uint128_t)prg.getRandom64();
	v = v >> 1;
	return ZpMersenne127Element(v);
}

