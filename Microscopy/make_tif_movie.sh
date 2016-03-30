#!/bin/bash

# Shell script to create a mp4 movie from processed images 
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Last Modified: March 02, 2016

# Movie framerate

FRAMESPERSEC="10"

# Convert tif to png and modify file name

for f in *.tif; do convert "$f" "$(basename "$f" .tif).png"; done
rename 's/_ch01//' *.png
rename 's/-60-40//' *.png
for f in *.png; do mv "$f" $(echo $f | tail -c 9); done

# Figure out filename

CNT="0"
NAME=""
for f in *.png; do
	if [ $CNT = "1" ]; then
		NAME=$(printf '%s' $(basename $f .png | sed 's/[0-9]*//g'))
	fi
	((CNT++))
done

LENGTH=$((CNT/FRAMESPERSEC))
echo "Creating Movie - Frames: $CNT Estimated Length: $LENGTH seconds"

# Encode movie

avconv -framerate $FRAMESPERSEC -i "$NAME%03d.png" -s:v 1280x720 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p "$NAME.mp4"
