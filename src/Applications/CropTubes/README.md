TubeTK CropTubes Application
=============================================

#### Overview:

Given a 3D Box and a TRE file containing a set of tubes, this module
returns tubes that are contained in the box.

#### Usage:

```
   CropTubes  [--boxCorner <std::vector<double>>]
              [--boxSize <std::vector<double>>]
              [--xml] [--echo] [--CropTubes]
              [--][--version] [-h]
              <std::string> <std::string>


Where:

   -p, --boxCorner <std::vector<double>>
     (required) 3D vector corresponding to the position of the bottom
     left-hand corner of the bounding box.

   -s, --boxSize <std::vector<double>>
     (required) 3D vector corresponding to the box dimensions
     along x,y and z axis.

   --xml
     Produce xml description of command line arguments (default: 0)

   --echo
     Echo the command line arguments (default: 0)

   -c, --CropTubes
     Specify if tubes should be cropped according to the box (default: 0)

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input TRE file containing tubes

   <std::string>
     (required)  Output TRE file containing only tubes that goes
     through the box


   Description: Return a set of tubes contained in a bounding box

   Author(s): Lucas Gandel (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.

```

#### Algorithm:

1. **Get the input tubes:**
 We read through the input TRE file in order to get each position of each
 tube.

2. **Determine if a tube belongs to the box:**
 For each position of the current tube, we check if the current point is
 located in the bounding box:

 * **If tubes should NOT be cropped** : We save the entire tube in a
 target group of tubes as long as one point is contained in the 3D box.

 * **If tubes should be cropped** : We read through all positions of
 the current tube. When a point that belongs to the box is found,
 its position is saved in the point list of a new target tube.
 If the next point still belongs to the box, it is saved in the same
 current target tube as the first one. On the contrary, if it does not
 belong to the box, the curent target tube is saved in the target group
 and we keep on checking the other positions. This enables to crop tubes
 that go in and out of the box.

3. **Write Output TRE file:**
 We write the target group of tubes in a TRE file.

----
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by
[Kitware, Inc.](http://www.kitware.com) and licensed under the
[Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*