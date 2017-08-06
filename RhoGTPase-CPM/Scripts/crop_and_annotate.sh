#! /bin/bash

for image in *.png; do

  # Crop image
  convert -crop 560X560+230+34 $image $image

  # Resize
  convert -resize 50% $image $image

  # Print frame number
  img_num=$(basename $image | sed -e 's/.png//' | tail -c 5)
  echo "Processing Image: $img_num"

  # Annotate
  convert $image -background Blue -fill White -pointsize 18\
  label:"MCS $img_num" -gravity Center -append $image

done
