/***** do_butterflies_neon_precalculated_twiddles.S *****/

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@ Do all the Butterflies for a stage using the neon unit @@@@@@@@@@@@@@@@@@@@@@@@@@ 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 	.syntax unified
 	.arch armv7-a
	.fpu neon
 	.section .text

	.align	2
	.global	do_butterflies_neon_precalculated_twiddles
	.thumb
	.thumb_func
	.type	do_butterflies_neon_precalculated_twiddles, %function
do_butterflies_neon_precalculated_twiddles:
	@ Perform all butterflies for one stage of the FFT.
	@ Inputs:	r0 = &fftTempBuffer				(inputs)
	@			r1 = &fftTempBufferInterleaved	(outputs)
	@			r2 = &TwiddleArray
	@
	@ Useful Registers:
	@			r3 = N (gFFTSize)
	
	vpush {q4, q5, q6, q7}
	
	@ Load the FFT Size.
	movw	r3, #:lower16:gFFTSize			@ Move the bottom 16 bits into r8, zero-ing the top 16 bits
	movt	r3, #:upper16:gFFTSize			@ Move the top 16 bits into r8
	ldr		r3, [r3, #0]
	asr   r3, r3, #1                  @ Do N/2 butterflies at each stage

	@ Perform a butterfly operation on two Complex inputs (a, b) and a Complex Twiddle factor (W) to give two Complex outputs (c, d).
	@ Inputs are assumed contiguous in memory so they can be accessed by offsets. Same with outputs.

	@ NEON Registers:
	@   q0 = d0, d1 = [W1.re, W2.re, W3.re, W4.re]
	@   q1 = d2, d3 = [W1.im, W2.im, W3.im, W4.im]
	
	@		q3 = d6,d7 = [4*a.re]
	@		q4 = d8,d9 = [4*a.im]
	@		q5 = d10,d11 = [4*b.re]
	@		q6 = d12,d13 = [4*b.im]
	@		q7 = d14,d15 = [4*P]
	@		q8 = d16,d17 = [4*Q]
	@		q9-12 = d18-d25 = results [c.re, c.im, d.re, d.im] to be stored.
	
do_butterfly_for_four_inputs:
	@ zero out the intermediate accumulators (d14-d17)
	veor q7, q7, q7
	veor q8, q8, q8

	@ Load the Twiddle factors in.
	vld2.32 {d0, d2}, [r2]!
	vld2.32 {d1, d3}, [r2]!
  
	@ Interleave load 4 inputs a and b
	vld4.32	{d6,d8,d10,d12}, [r0]!		@ '!' increases r0 by the amount of bytes read, in this case 2inputs*4floats*4bytes = 32 bytes.
	vld4.32 {d7,d9,d11,d13}, [r0]!
	
	@ Compute Intermediate results P and Q.
	vmla.f32	q7, q5, q0	   @ q7 = P = (W.re*b.re)
	vmls.f32	q7, q6, q1	   @					  - (W.im*b.im)
	vmla.f32	q8, q6, q0	   @ q8 = Q = (W.re*b.im)
	vmla.f32	q8, q5, q1	   @					  + (W.im*b.re)
	
	@ Compute Results
	vadd.f32	q9, q3, q7		@ q9 = c.re = a.re + P
	vadd.f32	q10, q4, q8		@ q10 = c.im = a.im + Q
	vsub.f32	q11, q3, q7		@ q11 = d.re = a.re - P
	vsub.f32	q12, q4, q8		@ q12 = d.im = a.im - Q
	
	@ Store the results into the results array.
	vst4.32 {d18, d20, d22, d24}, [r1]!
	vst4.32 {d19, d21, d23, d25}, [r1]!

	@add r7, r7, #4				@ Number of butterflies done for this type += 4.

	@Number of butterflies left to do -= 4
	subs r3, r3, #4
	@If we've not hit 0 yet, do 4 more!
	bne do_butterfly_for_four_inputs
	
done_all_butterfly_types:
	vpop {q4, q5, q6, q7}
	bx lr

