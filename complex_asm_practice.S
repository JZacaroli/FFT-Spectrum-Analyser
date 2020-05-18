@/***** butterfly_asm.S *****/

 	.syntax unified
 	.arch armv7-a
	.fpu vfpv3-d16
 	.section .text

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@ Complex Add @@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

@ void complexAdd_asm(Complex* a, Complex* b, Complex* c)
@
@ Registers:
@ 	r0: a*
@	r1: b*
@	r2: c*

	.align	2
	.global	complexAdd_asm
	.thumb
	.thumb_func
	.type	complexAdd_asm, %function
complexAdd_asm:
	@vpush {s4, r5, r6, r7}				@ Save registers before beginning
	
	vldr		s3, [r0, #0]		@ a.re
	vldr		s4,	[r0, #4]		@ a.im
	vldr		s5,	[r1, #0]		@ b.re
	vldr		s6,	[r1, #4]		@ b.im
	
	vadd.f32 s3, s3, s5				@ c.re = a.re + b.re
	vadd.f32 s4, s4, s6				@ c.im = a.im + b.im
	
	@ store the real part of the number into the address held by r2
	vstr s3, [r2, #0]
	@ and store the imaginary part of the number into the address held by r2 + #4
	vstr s4, [r2, #4]
	
	bx lr


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@ Complex Subtract @@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	.align	2
	.global	complexSubtract_asm
	.thumb
	.thumb_func
	.type	complexSubtract_asm, %function
complexSubtract_asm:
	
	vldr		s3, [r0, #0]		@ a.re
	vldr		s4,	[r0, #4]		@ a.im
	vldr		s5,	[r1, #0]		@ b.re
	vldr		s6,	[r1, #4]		@ b.im
	
	vsub.f32 s3, s3, s5				@ c.re = a.re - b.re
	vsub.f32 s4, s4, s6				@ c.im = a.im - b.im
	
	@ store the real and imaginary parts of the number in adjacent memory registers
	vstr s3, [r2, #0]
	vstr s4, [r2, #4]
	
	bx lr


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@ Complex Multiply @@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	.align	2
	.global	complexMultiply_asm
	.thumb
	.thumb_func
	.type	complexMultiply_asm, %function
complexMultiply_asm:
	
	vldr		s3, [r0, #0]		@ a.re
	vldr		s4,	[r0, #4]		@ a.im
	vldr		s5,	[r1, #0]		@ b.re
	vldr		s6,	[r1, #4]		@ b.im
	
	vmul.f32 s7, s3, s5				@ c.re = (a.re*b.re)
	vmls.f32 s7, s4, s6				@					 - (a.im*b.im)
									@
	vmul.f32 s8, s4, s5				@ c.im = (a.im*b.re) 
	vmla.f32 s8, s3, s6				@					 + (a.re*b.im)
	
	@ store the real and imaginary parts of the number in adjacent memory registers
	vstr s7, [r2, #0]
	vstr s8, [r2, #4]
	
	bx lr


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@ Butterfly @@@@@@@@@@@@@@@@@@@@@@@@@@ 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	.align	2
	.global	butterfly_asm
	.thumb
	.thumb_func
	.type	butterfly_asm, %function
butterfly_asm:
	@ Performs a butterfly operation on two Complex inputs (a, b) and a Complex Twiddle factor (W) to give two Complex outputs (c, d).
	@ Inputs are assumed contiguous in memory so they can be accessed by offsets. Same with outputs.
	@ r0: address of a.re
	@ r1: address of c.re
	@ r2: address of W.re
	@ s3-s6: Complex a,b. 
	@ s7,s8: Complex W
	@ s9,s10: Intermediate results P and Q
	@ s11-s14: Results to be stored into addresses pointed to by r1.
	
	@ Load inputs a and b
	vldr	s3, [r0, #0]		@ s3 = a.re
	vldr	s4,	[r0, #4]		@ s4 = a.im
	vldr	s5,	[r0, #8]		@ s5 = b.re
	vldr	s6,	[r0, #12]		@ s6 = b.im
	
	@ Load Twiddle factor W
	vldr	s7, [r2, #0]		@ s7 = W.re
	vldr	s8, [r2, #4]		@ s8 = W.im
	
	@ Compute Intermediate results P and Q.
	vmla.f32	s9, s7, s5		@ s9 = P = (W.re*b.re)
	vmls.f32	s9, s8, s6		@				      - (W.im*bim)
	
	vmla.f32	s10, s7, s6		@ s10 = Q = (W.re*b.im)
	vmla.f32	s10, s8, s5		@					   + (W.im*bre)
	
	@ Compute Results
	vadd.f32	s11, s3, s9		@ s3 = c.re = a.re + P
	vadd.f32	s12, s4, s10	@ s4 = c.im = a.im + Q
	vsub.f32	s13, s3, s9		@ s5 = d.re = a.re - P
	vsub.f32	s14, s4, s10	@ s6 = d.im = a.im - Q
	
	@ Store the results into the results vector.
	vstr s11, [r1, #0]
	vstr s12, [r1, #4]
	vstr s13, [r1, #8]
	vstr s14, [r1, #12]
	
	bx lr

