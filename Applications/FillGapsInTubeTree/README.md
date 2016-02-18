TubeTK FillGapsInTubeTree Application
=============================================

#### Overview:

Given an input tree, this application will interpolate the gap between parent and children in the tree.

#### Usage:

```
   FillGapsInTubeTree  [--method <std::string>]
              [--xml] [--echo]
              [--][--version] [-h]
              <std::string> <std::string>


Where:

   -method <Std::string>
     (required) The method use for interpolating the path to fill the gap in tube-tree.

   --xml
     Produce xml description of command line arguments (default: 0)

   --echo
     Echo the command line arguments (default: 0)

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input tre file containg tube structure in the form of a tree.

   <std::string>
     (required)  Output tre file


   Description: Return tube tree, where there is no gap between the parent and child.

   Author(s): Sumedha Singla (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.

```

#### Algorithm:

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