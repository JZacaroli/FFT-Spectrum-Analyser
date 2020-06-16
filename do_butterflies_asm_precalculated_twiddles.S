/***** do_butterflies_asm_precalculated_twiddles.S *****/
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@ Do all the Butterflies for an FFT stage @@@@@@@@@@@@@@@@@@@@@@@@@@ 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 	.syntax unified
 	.arch armv7-a
	.fpu vfpv3-d16
 	.section .text

	.align	2
	.global	do_butterflies_asm_precalculated_twiddles
	.thumb
	.thumb_func
	.type	do_butterflies_asm_precalculated_twiddles, %function
do_butterflies_asm_precalculated_twiddles:
	@ Perform all butterflies for one stage of the FFT.
	@ Inputs:	r0 = &fftTempBuffer				(inputs)
	@			r1 = &fftTempBufferInterleaved	(outputs)
	@			r2 = &twiddleArray
	@
	@ Useful Registers:
	@			r3 = N (gFFTSize)
	
	@ Load the FFT Size.
	movw	r3, #:lower16:gFFTSize			@ Move the bottom 16 bits into r3, zero-ing the top 16 bits
	movt	r3, #:upper16:gFFTSize			@ Move the top 16 bits into r3
	ldr		r3, [r3, #0]
	@ Number of butterflies to do: r3 = N/2
	asr		r3, r3, #1

	vpush {s4, s5, s6, s7, s8, s9, s10, s11}			@ Save registers before beginning.

	@ Perform a butterfly operation on two Complex inputs (a, b) and a Complex Twiddle factor (W) to give two Complex outputs (c, d).
	@ Inputs are assumed contiguous in memory so they can be accessed by offsets. Same with outputs.
	@		s0 = W.re
	@		s1 = W.im
	@		s2 = a.re
	@		s3 = a.im
	@		s4 = b.re
	@		s5 = b.im
	@		s6 = P
	@		s7 = Q
	@		s8-s11 = results to be stored.
	
do_single_butterfly:
	@ zero out the intermediate accumulators
	vmov s6, r10
	vmov s7, r10
	
	@ Load inputs a and b
	vldr	s2, [r0, #0]		@ s2 = a.re
	vldr	s3,	[r0, #4]		@ s3 = a.im
	vldr	s4,	[r0, #8]		@ s4 = b.re
	vldr	s5,	[r0, #12]		@ s5 = b.im
	
	@Load the twiddle factor
	vldr	s0, [r2, #0]
	vldr	s1, [r2, #4]
	
	@ Compute Intermediate results P and Q.
	vmla.f32	s6, s0, s4		@ s6 = P = (W.re*b.re)
	vmls.f32	s6, s1, s5		@				      - (W.im*bim)
	
	vmla.f32	s7, s0, s5		@ s7 = Q = (W.re*b.im)
	vmla.f32	s7, s1, s4		@					   + (W.im*bre)
	
	@ Compute Results
	vadd.f32	s8, s2, s6		@ s8 = c.re = a.re + P
	vadd.f32	s9, s3, s7		@ s9 = c.im = a.im + Q
	vsub.f32	s10, s2, s6		@ s10 = d.re = a.re - P
	vsub.f32	s11, s3, s7		@ s11 = d.im = a.im - Q
	
	@ Store the results into the results array.
	vstr s8, [r1, #0]
	vstr s9, [r1, #4]
	vstr s10, [r1, #8]
	vstr s11, [r1, #12]
	
	subs r3, r3, #1				@ Number of butterflies left to do -= 1
	add r0, r0, #16				@ Move where we're indexing into each fft buffer forward by 16 bytes. (4 4-byte floats)
	add r1, r1, #16
	add r2, r2, #8				@ Move where we're indexing into the twiddle factor array forward by 8 bytes (2 4-byte floats)
	
	@ If the number of butterflies left to do is greater than 0, go to do_single_butterfly
	bne do_single_butterfly
	
done_all_butterflies:
	vpop {s4, s5, s6, s7, s8, s9, s10, s11}		@ Reset registers at end.
	bx lr

