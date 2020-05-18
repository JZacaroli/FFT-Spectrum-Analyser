/***** complex_jz.h *****/
/**
 * Includes all the necessary methods to manipulate complex numbers for the FFT.
 * 
 * Would probably be nicer as an object with overloaded operators. Oh well.
 */

/**
 * Complex Number struct
 */
struct Complex {
	float re=1, im=0;
};

/* Assembly-language Complex add function */
extern "C" {
	void complexAdd_asm(Complex* a, Complex* b, Complex* c);
	void complexSubtract_asm(Complex* a, Complex* b, Complex* c);
	void complexMultiply_asm(Complex* a, Complex* b, Complex* c);
	void butterfly_asm(Complex* inputs, Complex* outputs, Complex* twiddle);
}

/**
 * complexMultiply - Multiply two complex numbers
 */
Complex complexMultiply(Complex a, Complex b) {
	Complex c;
	c.re = (a.re*b.re) - (a.im*b.im);
	c.im = (a.re*b.im) + (a.im*b.re);
	return c;
}

/**
 * complexAdd - Add two complex numbers
 */
Complex complexAdd(Complex a, Complex b) {
	Complex c;
	c.re = a.re + b.re;
	c.im = a.im + b.im;
	return c;
}

/**
 * complexSubtract - Subtract a complex number (b) from another (a)
 */
Complex complexSubtract(Complex a, Complex b) {
	Complex c;
	c.re = a.re-b.re;
	c.im = a.im-b.im;
	return c;
}

/**
 * calculateNegativeComplex - Finds the negative of a complex number struct
 */
Complex calculateNegativeComplex(Complex a) {
	Complex b;
	b.re = -a.re;
	b.im = -a.im;
	return b;
}