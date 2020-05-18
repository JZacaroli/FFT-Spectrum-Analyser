/***** do_butterflies_asm.S *****/

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@ Do all the Butterflies for an FFT stage @@@@@@@@@@@@@@@@@@@@@@@@@@ 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 	.syntax unified
 	.arch armv7-a
	.fpu vfpv3-d16
 	.section .text

	.align	2
	.global	do_butterflies_asm
	.thumb
	.thumb_func
	.type	do_butterflies_asm, %function
do_butterflies_asm:
	@ Perform all butterflies for one stage of the FFT.
	@ Inputs:	r0 = &fftTempBuffer				(inputs)
	@			r1 = &fftTempBufferInterleaved	(outputs)
	@			r2 = &twiddleRotationFactor		(used to rotate the twiddle at each stage)
	@			r3 = (int) stage
	@
	@ Useful Registers:
	@			r4 = numOfTypesOfButterfliesToDo (int)
	@			r5 = numOfTypesOfButterfliesDone (int)
	@			r6 = numOfButterfliesToDoForType (int)
	@			r7 = numOfButterfliesDoneForType (int)
	@			r8 = N (gFFTSize)
	@
	
	push {r4, r5, r6, r7, r8, r9, r10}		@ Save registers before beginning
	
	@ Load the FFT Size.
	movw	r8, #:lower16:gFFTSize			@ Move the bottom 16 bits into r8, zero-ing the top 16 bits
	movt	r8, #:upper16:gFFTSize			@ Move the top 16 bits into r8
	ldr		r8, [r8, #0]

	@ There's 2^stage different types of butterfly in each stage
	mov		r4, #1
	LSL		r4, r4, r3			@ left shift #1 by stage. = 2^stage.
	mov r5, #0					@ numOfTypesOfButterfliesDone=0
	
	
	@ There's N/(2^(stage+1)) of each butterfly type.
	mov		r6, r8
	ASR 	r6, r6, #1			@ Divide N by 2
	ASR		r6, r6, r3			@ And then divide by 2^stage.
	
	@ Initial Twiddle Factor = (1 + j0)
	mov		r9, #1
	mov 	r10, #0
	vmov	s0, s1, r9, r10
	vcvt.f32.s32 s0, s0
	vcvt.f32.s32 s1, s1

	@ Load the twiddle rotation factor onto the s registers
	vpush {s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15}			@ Save registers before beginning. TODO: DO WE NEED TO DO THIS?!?!?
	vldr	s12, [r2, #0]		@ W_rot.re
	vldr	s13, [r2, #4]		@ W_rot.im

do_butterfly_type:
	mov r7, #0					@ Number of butterflies done for this type = 0.

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
	@		s12 = W_rot.re   Twiddle rotation factor W_rot. Used after each butterfly type.
	@		s13 = W_rot.im
	@		s14,s15 = intermediate result registers for rotating twiddle factor.
do_single_butterfly:
	@ zero out the intermediate accumulators
	vmov s6, r10
	vmov s7, r10
	
	@ Load inputs a and b
	vldr	s2, [r0, #0]		@ s2 = a.re
	vldr	s3,	[r0, #4]		@ s3 = a.im
	vldr	s4,	[r0, #8]		@ s4 = b.re
	vldr	s5,	[r0, #12]		@ s5 = b.im
	
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
	
	add r7, r7, #1				@ Number of butterflies done for this type += 1.
	add r0, r0, #16				@ Move where we're indexing into each fft buffer forward by 16 bytes. (4 4-byte floats)
	add r1, r1, #16
	
	@ If we've not done all the butterflies for this type, go to do_single_butterfly
	cmp r6, r7
	bne do_single_butterfly
	
	@ We've done all butterflies for this type
	add r5, r5, #1				@numOfTypesOfButterfliesDone+=1
	
	@ If we've done all the butterfly types, finish!
	cmp r4, r5
	beq done_all_butterfly_types
	
	@ If we've not done all the butterfly types, rotate the twiddle factor (W) by the twiddleRotationFactor (W_rot) (i.e. W*=W_rot)
	vmul.f32	s14, s0, s12			@ s14 = W.re = (W.re*W_rot.re)
	vmls.f32	s14, s1, s13			@							  - (W.im*W_rot.im)
	vmul.f32	s15, s1, s12			@ s15 = W.im = (W.im*W_rot.re) 
	vmla.f32	s15, s0, s13			@							  + (W.re*W_rot.im)
	vmov.f32	s0, s14					@ s0 = W.re
	vmov.f32	s1, s15					@ s1 = W.im
	
	@ Then go to do the next butterfly type.
	b do_butterfly_type
	
done_all_butterfly_types:
	vpop {s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15}
	pop {r4, r5, r6, r7, r8, r9, r10}		@ Reset registers at end.
	bx lr

