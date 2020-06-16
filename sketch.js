/**
 * Messages to be sent back and forth between GUI and render.cpp:
 *			Buffer 0: Whether or not filtering is happening [0] / [1]
 *			Buffer 1: What the filtering gain values are [gainLow, gainMid, gainHigh]
 *			Buffer 2: Receiving the FFT data on this buffer. [k0, k1, ...., kN/2]
 */

var guiSketch = new p5(function( sketch ) {
	
	class GainSlider {
	  constructor(x, y, w, h, r, safeHeight) {
	    this.dragging = false; // Is the object being dragged?
	    this.rollover = false; // Is the mouse over the ellipse?
	    this.start_x = x;
	    this.start_y = y;
	    this.safeHeight = safeHeight-h;
	    this.x = x;
	    this.y = y;
	    this.w = w;
	    this.h = h;
	    this.r = r;
	    this.offsetX = 0;
	    this.offsetY = 0;
	    this.gain = 1;
	  }
	
	  over() {
	    // Is mouse over object
	    if (sketch.mouseX > this.x && sketch.mouseX < this.x + this.w && sketch.mouseY > this.y && sketch.mouseY < this.y + this.h) {
	      this.rollover = true;
	    } else {
	      this.rollover = false;
	    }
	  }
	
	  update() {
	    // Adjust location if being dragged
	    if (this.dragging) {
	    	if (sketch.mouseY < this.start_y + this.h) {
	    		this.y = this.start_y;
	    		this.gain = 1;
	    	} else if (sketch.mouseY > this.start_y+this.safeHeight) {
	    		this.y = this.start_y + this.safeHeight;
	    		this.gain = 0;
	    	} else {
	    		this.y = sketch.mouseY + this.offsetY;
	    		// this.gain = 0.5;
	    	}
	    	this.gain = 1 - (this.y-this.start_y)/this.safeHeight;
	    }
	  }
	
	  show() {
	    sketch.stroke(0);
	    // Different fill based on state
	    if (this.dragging) {
	      sketch.fill(50);
	    } else if (this.rollover) {
	      sketch.fill(100);
	    } else {
	      sketch.fill(175, 200);
	    }
	    sketch.rect(this.x, this.y, this.w, this.h, this.r);
	  }
	
	  pressed() {
	    // Did I click on the rectangle?
	    if (sketch.mouseX > this.x && sketch.mouseX < this.x + this.w && sketch.mouseY > this.y && sketch.mouseY < this.y + this.h) {
	      this.dragging = true;
	      // If so, keep track of relative location of click to corner of rectangle
	      this.offsetX = this.x - sketch.mouseX;
	      this.offsetY = this.y - sketch.mouseY;
	    }
	  }
	
	  released() {
	    // Quit dragging
	    this.dragging = false;
	  }
	}
    let canvas_dimensions = [sketch.windowWidth, sketch.windowHeight];
    
    let sampleRate = 44100;
    var frameCount = 0;
    var sendGainBufferInterval = 1;
    
    //Spectrum Display variables
	let FftDisplayXBegin = 50;
	let FftDisplayYBegin = 40;
	let FftDisplayWidth = 600;
	let FftDisplayHeight = 300;
	let FftDisplayInsideMargin = 0;
	let FftDisplayWidthToPlayWith = FftDisplayWidth - 2*FftDisplayInsideMargin;
	
	// EQ button
	var FilteringOn = false;
	let filteringButtonXBegin = 125;
	let filteringButtonWidth = 50;
	let filteringButtonHeight = 50;
	let filteringButtonYBegin = FftDisplayYBegin + FftDisplayHeight + 70;
	
	// EQ Sliders
	let marginSliders = 70;
	let sliderWidth = 50;
	let sliderHeight = 150;
	//Low
	let sliderLowXBegin = filteringButtonXBegin+filteringButtonWidth+marginSliders;
	let sliderLowYBegin = filteringButtonYBegin;
	//Mid
	let sliderMidXBegin = sliderLowXBegin + sliderWidth + marginSliders;
	let sliderMidYBegin = filteringButtonYBegin;
	//High
	let sliderHighXBegin = sliderMidXBegin + sliderWidth + marginSliders;
	let sliderHighYBegin = filteringButtonYBegin;
	
    sketch.setup = function() {
        sketch.createCanvas(canvas_dimensions[0], canvas_dimensions[1]);
        sketch.textFont('Courier New');
        
        lowGainSlider = new GainSlider(sliderLowXBegin, sliderLowYBegin, sliderWidth, 15, 2, sliderHeight);
        midGainSlider = new GainSlider(sliderMidXBegin, sliderMidYBegin, sliderWidth, 15, 2, sliderHeight);
        highGainSlider = new GainSlider(sliderHighXBegin, sliderHighYBegin, sliderWidth, 15, 2, sliderHeight);
    };
    
    sketch.mousePressed = function() {
    	if (sketch.mouseX > filteringButtonXBegin && sketch.mouseX < filteringButtonXBegin+filteringButtonWidth && sketch.mouseY > filteringButtonYBegin && sketch.mouseY < filteringButtonYBegin+filteringButtonHeight) {
    		FilteringOn = !FilteringOn;
    		if (FilteringOn) {
    			Bela.data.sendBuffer(0, 'int', 1);
    		} else {
    			Bela.data.sendBuffer(0, 'int', 0);
    		}
    	}
    	lowGainSlider.pressed();
    	midGainSlider.pressed();
    	highGainSlider.pressed();
    	
    }
    
    sketch.mouseReleased = function() {
    	lowGainSlider.released();
    	midGainSlider.released();
    	highGainSlider.released();
    }
	
    sketch.draw = function() {
    	frameCount++;
    	
    	if (FilteringOn) {
			if (frameCount%sendGainBufferInterval === 0) {
				Bela.data.sendBuffer(1, 'float', [lowGainSlider.gain, midGainSlider.gain, highGainSlider.gain]);
			}
		}
    	
        sketch.background(250);
		
		//------------------------//
        // Draw Spectrum Analyzer //
        //------------------------//
        
        sketch.fill(51, 153, 255);
        sketch.rect(FftDisplayXBegin, FftDisplayYBegin, FftDisplayWidth, FftDisplayHeight);
        sketch.stroke(0);
        sketch.fill(0);
        sketch.textAlign(sketch.CENTER);
        sketch.textSize(30);
		sketch.text('Spectrum Display', FftDisplayXBegin+FftDisplayWidth/2, 25);
        sketch.textSize(15);
        sketch.text('20Hz', FftDisplayXBegin, FftDisplayYBegin+FftDisplayHeight+15);
        sketch.text('200Hz', FftDisplayXBegin + FftDisplayWidth/3, FftDisplayYBegin+FftDisplayHeight+15);
        sketch.text('2000Hz', FftDisplayXBegin + 2*FftDisplayWidth/3, FftDisplayYBegin+FftDisplayHeight+15);
        sketch.text('20000Hz', FftDisplayXBegin + FftDisplayWidth, FftDisplayYBegin+FftDisplayHeight+15);
        
        if (FilteringOn) {
        	sketch.fill(255);
        } else {
        	sketch.fill(150);
        }
        sketch.rect(filteringButtonXBegin-20, filteringButtonYBegin-40, 4*(sliderWidth+marginSliders), 1.4*sliderHeight, 10, 10, 10, 10);
        
        //-------------------------//
        // Filtering on/off button //
        //-------------------------//
        //Filtering button/text colours all defined by this.
        if (FilteringOn) {
        	sketch.fill(50, 100, 230);
        } else {
        	sketch.fill(100, 100, 100);
        }
        sketch.textSize(16);
        sketch.text('EQ', filteringButtonXBegin+filteringButtonWidth/2, filteringButtonYBegin-10);
        sketch.rect(filteringButtonXBegin, filteringButtonYBegin, filteringButtonWidth, filteringButtonHeight, 10, 10, 10, 10);
        
        //------------------//
        // Draw the sliders //
        //------------------//
        sketch.text('Lows', sliderLowXBegin+sliderWidth/2, filteringButtonYBegin-10);
        sketch.text('Mids', sliderMidXBegin+sliderWidth/2, filteringButtonYBegin-10);
        sketch.text('Highs', sliderHighXBegin+sliderWidth/2, filteringButtonYBegin-10);
        sketch.rect(sliderLowXBegin+(sliderWidth/2)-2, sliderLowYBegin, 4, sliderHeight, 2, 2, 2, 2);
        sketch.rect(sliderMidXBegin+(sliderWidth/2)-2, sliderMidYBegin, 4, sliderHeight, 2, 2, 2, 2);
        sketch.rect(sliderHighXBegin+(sliderWidth/2)-2, sliderHighYBegin, 4, sliderHeight, 2, 2, 2, 2);
        lowGainSlider.over();
        lowGainSlider.update();
        lowGainSlider.show();
        midGainSlider.over();
        midGainSlider.update();
        midGainSlider.show();
        highGainSlider.over();
        highGainSlider.update();
        highGainSlider.show();
        
        //------------------------------//
        // Draw lines over the fft plot //
        //------------------------------//
        sketch.strokeWeight(1);
        sketch.stroke(150);
        for (var i=0; i<6; i++) {
        	sketch.line(FftDisplayXBegin, FftDisplayYBegin+(i*FftDisplayHeight/5), FftDisplayXBegin+FftDisplayWidth, FftDisplayYBegin+(i*FftDisplayHeight/5));
        }
        for (var i=0; i<4; i++) {
        	sketch.line(FftDisplayXBegin+(i*FftDisplayWidth/3), FftDisplayYBegin, FftDisplayXBegin+(i*FftDisplayWidth/3), FftDisplayYBegin+FftDisplayHeight);
        }
        
        
        //--------------------------//
        // Get FFT Data and plot it //
        //--------------------------//
        let fftData = Bela.data.buffers[2];
        fftSize = fftData.length; //This is only half the length of the actual FFT because I'm only sending the first half of the spectrum.
        
        //Only do frequencies between 20Hz-20kHz
        var start_k=0;
        while (start_k*sampleRate/(2*fftSize) < 20) {
        	start_k++;
        }
        var end_k=fftData.length-1;
        while (end_k*sampleRate/(2*fftSize) > 20000) {
        	end_k--;
        }
        
        sketch.strokeWeight(1);
    	sketch.stroke(250);
    	var prevLineX = 0;
    	var prevLineHeight = -1;
    	for (var i = start_k; i<end_k; i++) {
    		//let lineXBegin = i*FftDisplayWidthToPlayWith/(fftSize-1);
	    	var lineFrequency = i*sampleRate/(2*fftSize);
    		if (lineFrequency>20 && lineFrequency<20000)
    		{
	    		//let lineX = (Math.log(i-4-start_k)/Math.log(10)) * FftDisplayWidthToPlayWith / (Math.log(end_k-start_k)/Math.log(10));
	    		let lineX = (Math.log(lineFrequency/20)/Math.log(10)) * FftDisplayWidthToPlayWith/3;
	    		var lineHeight = -Math.log(0.00001+fftData[i])/Math.log(10); //Between 0 and 5
	    		lineHeight = (FftDisplayHeight*lineHeight/5); //Between 0 and display height
	    		if (prevLineHeight == -1) {
	    			prevLineHeight = lineHeight;
	    		}
	    		sketch.line(FftDisplayXBegin+FftDisplayInsideMargin + prevLineX,
	    					FftDisplayYBegin+prevLineHeight,
	    					FftDisplayXBegin+FftDisplayInsideMargin + lineX,
	    					FftDisplayYBegin+lineHeight);
	    		prevLineX = lineX;
	    		prevLineHeight = lineHeight;
    		}
    	}
    }
    
}, 'gui');

