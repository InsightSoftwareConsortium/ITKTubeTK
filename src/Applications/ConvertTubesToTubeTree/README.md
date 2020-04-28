TubeTK ConvertTubesToTubeTree Application
=============================================

#### Overview:

Given a TRE file containing a disjoint set of tubes, this module tries to
 compute the connectivity between tubes.

#### Usage:

```
   ConvertTubesToTubeTree  [--returnparameterfile <std::string>]
                           [--processinformationaddress <std::string>]
                           [--xml] [--echo] [--removeOrphanTubes]
                           [--maxContinuityAngleError <double>]
                           [--maxTubeDistanceToRadiusRatio <double>]
                           [--rootTubeIdList <std::vector<int>>] [--]
                           [--version] [-h] <std::string> <std::string>


Where:

   --returnparameterfile <std::string>
     Filename in which to write simple return parameters (int, float,
     int-vector, etc.) as opposed to bulk return parameters (image,
     geometry, transform, measurement, table).

   --processinformationaddress <std::string>
     Address of a structure to store process information (progress, abort,
     etc.). (default: 0)

   --xml
     Produce xml description of command line arguments (default: 0)

   --echo
     Echo the command line arguments (default: 0)

   --removeOrphanTubes
     Specify if orphan tubes should be removed (default: 0)

   --maxContinuityAngleError <double>
     Maximum continuity angle error (default: 180)

   --maxTubeDistanceToRadiusRatio <double>
     Maximum allowed distance to radius ratio to connect tubes (default: 2)

   --rootTubeIdList <std::vector<int>>
     List of root tube Ids

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input TRE file containing disjoint tubes

   <std::string>
     (required)  Output TRE file containing a tree of tubes


   Description: Constructs a tree from a set of disjoint tubes

   Author(s): Deepak Roy Chittajallu (Kitware), Matt McCormick (Kitware),
     Stephen Aylward (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.

```

#### Algorithm:

1. **Build tube adjancency graph:**
 Given a tre file with a set of tubes, we build a weighted directed tube
 adjacency graph wherein vertices correspond to tubes in the given tre file
 and edges connect adjacent tubes. We add an edge from tube A to tube B in the
 tube adjacency graph, if the following two criterion are satisfied for any of
 the two end points of tube B:

  * **Distance criterion**: distance of the end point of tube B to its closest
 point in tube A is less than *maxTubeDistanceToRadiusRatio* times the radius
 of the closest point in tube A.

  * **Continuity angle criterion**: The angle between the "vector ending at the
  end point of tube B and starting at its closest point in tube A" and the
  "vector joining the end point of tube B to its adjacent point in tube B" is
  less than *maxContinuityAngleError*.

 The weights of the edges are set equal to the minimum of the distance between
 either end points of tube B to tube A that satisfy the above two criterion.

2. **Prepare a list of root tubes:**
 Since the tube adjacency graph could potentially contain multiple connected
 components, we allow the user to choose the desired connected components by
 specifiying the *rootTubeIdList* parameter. If this is not provided, we
 consider each input tube as a potential root to retain all connected
 components in output.

3. **Build tube trees:**
 We pick root tubes in decreasing order of out-degree in the tube adjacency
 graph (in case of duplicates, pick in decreasing order of length) and run a
 minimum spanning tree style algorithm to obtain a tree of all tubes reachable
 from it in the tube adjacency graph. If the list of root tube ids contains
 multiple root tubes within the same connected component, then the root tube
 with largest out-degree with largest number of points will end up being the
 root for its connected component and all others will become its descendants in
 the generated tube tree.

----
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by
[Kitware, Inc.](http://www.kitware.com) and licensed under the
[Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*