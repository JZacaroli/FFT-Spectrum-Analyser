/***** do_butterflies_neon.S *****/

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@ Do all the Butterflies for a stage using the neon unit @@@@@@@@@@@@@@@@@@@@@@@@@@ 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 	.syntax unified
 	.arch armv7-a
	.fpu neon
 	.section .text

	.align	2
	.global	do_butterflies_neon
	.thumb
	.thumb_func
	.type	do_butterflies_neon, %function
do_butterflies_neon:
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
	@			r9 = 1
	@			r10 = 0
	
	
	push {r4, r5, r6, r7, r8, r9, r10}		@ Save registers before beginning
	vpush {q4, q5, q6, q7}
	
	
	@ Load the FFT Size.
	movw	r8, #:lower16:gFFTSize			@ Move the bottom 16 bits into r8, zero-ing the top 16 bits
	movt	r8, #:upper16:gFFTSize			@ Move the top 16 bits into r8
	ldr		r8, [r8, #0]

	@ There's r4 = 2^stage different types of butterfly in each stage
	mov		r4, #1
	LSL		r4, r4, r3			@ left shift #1 by stage. = 2^stage.
	mov r5, #0					@ numOfTypesOfButterfliesDone=0
	
	@ There's r6 = N/(2^(stage+1)) of each butterfly type.
	mov		r6, r8
	ASR 	r6, r6, #1			@ Divide N by 2
	ASR		r6, r6, r3			@ And then divide by 2^stage.
	
	@ Initial Twiddle Factor d0 = {s0,s1} = (1 + j0)
	mov		r9, #1
	mov 	r10, #0
	vmov	d0[0], r9
	vmov	d1[0], r10
	vcvt.f32.s32 d0, d0		@ Convert the int to a float representation
	vcvt.f32.s32 d1, d1		@ Convert int to float.

	@ Load the twiddle rotation factor at address=r2 onto the d registers d2, d3
	vld1.32 {d2}, [r2]
	vext.32 d3, d2, d2, #1		@ This seems an ugly way of doing this..
	
do_butterfly_type:
	mov r7, #0					@ Number of butterflies done for this type = 0.

	@ Perform a butterfly operation on two Complex inputs (a, b) and a Complex Twiddle factor (W) to give two Complex outputs (c, d).
	@ Inputs are assumed contiguous in memory so they can be accessed by offsets. Same with outputs.

	@ NEON Registers:
	@		d0 = W.re
	@		d1 = W.im
	@		d2 = W_rot.re
	@		d3 = W_rot.im
	@		d4,d5 = Intermediate results for calculating rotation of twiddles.
	
	@		q3 = d6,d7 = [4*a.re]
	@		q4 = d8,d9 = [4*a.im]
	@		q5 = d10,d11 = [4*b.re]
	@		q6 = d12,d13 = [4*b.im]
	@		q7 = d14,d15 = [4*P]
	@		q8 = d16,d17 = [4*Q]
	@		q9-12 = d18-d25 = results [c.re, c.im, d.re, d.im] to be stored.
	
do_single_butterfly_for_four_inputs:
	@ zero out the intermediate accumulators (d14-d17)
	veor q7, q7, q7
	veor q8, q8, q8
	
	@ Interleave load 4 inputs a and b that all go through the same butterfly.
	vld4.32	{d6,d8,d10,d12}, [r0]!		@ '!' increases r0 by the amount of bytes read, in this case 2inputs*4floats*4bytes = 32 bytes.
	vld4.32 {d7,d9,d11,d13}, [r0]!
	
	
	@ Compute Intermediate results P and Q.
	vmla.f32	q7, q5, d0[0]	@ q7 = P = (W.re*b.re)
	vmls.f32	q7, q6, d1[0]	@					  - (W.im*b.im)
	vmla.f32	q8, q6, d0[0]	@ q8 = Q = (W.re*b.im)
	vmla.f32	q8, q5, d1[0]	@					  + (W.im*b.re)
	
	@ Compute Results
	vadd.f32	q9, q3, q7		@ q9 = c.re = a.re + P
	vadd.f32	q10, q4, q8		@ q10 = c.im = a.im + Q
	vsub.f32	q11, q3, q7		@ q11 = d.re = a.re - P
	vsub.f32	q12, q4, q8		@ q12 = d.im = a.im - Q
	
	@ Store the results into the results array.
	vst4.32 {d18, d20, d22, d24}, [r1]!
	vst4.32 {d19, d21, d23, d25}, [r1]!
	
	add r7, r7, #4				@ Number of butterflies done for this type += 4.
	
	@ If we've not done all the butterflies for this type, go to do_single_butterfly
	cmp r6, r7
	bne do_single_butterfly_for_four_inputs
	
	@ We've done all butterflies for this type
	add r5, r5, #1				@numOfTypesOfButterfliesDone+=1
	
	
	@ If we've done all the butterfly types, finish!
	cmp r4, r5
	beq done_all_butterfly_types
	
	@ If we've not done all the butterfly types, rotate the twiddle factor (W) by the twiddleRotationFactor (W_rot) (i.e. W*=W_rot)
	vmul.f32	d4, d0, d2				@ d4 = W.re = (W.re*W_rot.re)
	vmls.f32	d4, d1, d3				@							  - (W.im*W_rot.im)
	vmul.f32	d5, d1, d2				@ d5 = W.im = (W.im*W_rot.re) 
	vmla.f32	d5, d0, d3				@							  + (W.re*W_rot.im)
										@ d0 and d1 are the rightful heirs to W.re and W.im
	vmov.f32	d0, d4					@ d0 = W.re
	vmov.f32	d1, d5					@ d1 = W.im
	
	@ Then go to do the next butterfly type.
	b do_butterfly_type
	
done_all_butterfly_types:
	vpop {q4, q5, q6, q7}
	pop {r4, r5, r6, r7, r8, r9, r10}		@ Reset registers at end.
	bx lr

