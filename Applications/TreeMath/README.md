TubeTK TreeMath Application
=============================================

#### Overview:

Given an input tree, this application will will perform different functions on tree like
interpolating the gap between parent and children in the tree.

#### Usage:

```
	 TreeMath
   Option infile is required but not defined
 Command tags:
   [ -w < filename > ]
      = Writes current tubes to the designated file.
        With: filename = Output filename
   [ -f < InterpolationMethod > ]
      = Connects the parent and child tube if they have a gap inbetween, by interpolating the path inbetween.
        With: InterpolationMethod = [S]traight Line, [L]Linear Interpolation, [C]urve Fitting, [M]inimal Path
 Command fields:
   < infile >
      = infile filename

   Description: Return tube tree, where there is no gap between the parent and child.

   Author(s): Sumedha Singla (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.

```

#### Algorithm for filling gap:

1. **Get the input tube tree:**
 The first step is to read the input tube-tree. The input tube-tree is a tre file, usually the output of
 application "ConvertTubesToTubeTree".

2. **Find the tubes with gap**
 For each parent tube we will iterate through all its children and find the nearest point
 on child tube from the parent's tubes end points.
 If the distance between these points is more than a threshold, there is a gap.

 * **Fill the gap**
 To fill the gap we will interpolate the gap between these two points using the method choosen by the user.

3. **Write Output TRE file**
 The final tube tree is writen as TRE file.

----
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by
[Kitware, Inc.](http://www.kitware.com) and licensed under the
[Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*