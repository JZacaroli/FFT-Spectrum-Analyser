/***** bitReversal.S *****/
@ Simple assembly function for bit reversal
@ Idea taken from Thomas Matthews' answer:
@  https://stackoverflow.com/questions/2602823/in-c-c-whats-the-simplest-way-to-reverse-the-order-of-bits-in-a-byte
@ The basic idea is that we're rotating one register's bits out from the bottom end
@ and into another register (that rotates the opposite way) from the bottom end.
@
@ Joe Zacaroli

@ Declare unified (ARM and Thumb) syntax
 	.syntax unified
 	.section .text

@ Functions pass the first four arguments on r0, r1, r2 and r3
@ Return values come back on r0.
@ r12 can be modified without restoring at the end
@ Don't touch r13-r15.
@ Other registers can be modified if restored at the end

@ Use the following lines before ANY function:
	.align 2						@ Align the function code to a 4-byte (2^2) word boundary.
	.arm							@ Use ARM instructions instead of Thumb.
	.type bitReversal, %function	@ Declare that this is a function
	.global bitReversal 			@ Set this function to be global (accessible from outside)
bitReversal:						@ Declare the function
	@ int to reverse is in r0
	@ number of bits to reverse by is in r1
	
	@ zero out the r2 register.
	mov r2, #0
	
rotateAnotherBit:
	RORS r0, r0, #1			@ Rotate r0 right by one bit and set the carry flag.
	LSL r2, r2, #1			@ Shift r2 left by 1 bit.
	addcs r2, r2, #1		@ Add 1 to r2 if carry is set.
	subs r1, r1, #1			@ Decrement the number of bits to reverse by 1 and compare to 0.
	bgt rotateAnotherBit	@ If num of bits to reverse is greater than 0, rotateAnotherBit!
	
	mov r0, r2				@ Move the second register to the return register.
	bx lr					@ Return (return value is in r0)
