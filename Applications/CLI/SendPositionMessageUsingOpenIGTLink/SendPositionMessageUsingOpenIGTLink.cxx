/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved. 

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// The following three should be used in every CLI application
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"
#include "tubeMessage.h"

// Includes specific to this CLI application
#include "igtlPositionMessage.h"
#include "igtlServerSocket.h"
#include "igtlOSUtil.h"

// Must include CLP before including Helper functions
#include "SendPositionMessageUsingOpenIGTLinkCLP.h"

void GenerateRandomPosition( double position[3] );
void GenerateRandomOrientation( double quaternion[4] );

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  int interval = (int) (1000.0 / frequency);
  std::cout << "Position information generation interval:\t" << interval << std::endl;

  igtl::PositionMessage::Pointer positionMsg;
  positionMsg = igtl::PositionMessage::New();
  positionMsg->SetDeviceName("Tracker");

  igtl::ServerSocket::Pointer serverSocket;
  serverSocket = igtl::ServerSocket::New();
  int serverCreationStatus = serverSocket->CreateServer(port);

  if (serverCreationStatus < 0)
    {
    std::cerr << "Cannot create a server socket." << std::endl;
    return EXIT_FAILURE;
    }

  igtl::Socket::Pointer socket;
  
  int time = 0;
  while (time < duration || duration < 0)
    {
    // Waiting for Connection
    std::cout << "Waiting for connection..." << std::endl;
    socket = serverSocket->WaitForConnection(1000);
    time++;
    
    if (socket.IsNotNull()) // if client connected
      {
      std::cout << "\t Got a connection with a client socket.." << std::endl;

      for (int i=0; i < 100; i++ )
        {
        double position[3];
        GenerateRandomPosition(position);

        double quaternion[4];
        GenerateRandomOrientation( quaternion );

        std::cout << "Position\t(" 
                  << position[0] << "," 
                  << position[1] << "," 
                  << position[2] << ")" << std::endl;
        std::cout << "Orientation\t" 
                  << quaternion[0] << "," 
                  << quaternion[1] << "," 
                  << quaternion[2] << ")" << std::endl;

        positionMsg->SetPosition( position[0], position[1], position[2] );
        positionMsg->SetQuaternion( quaternion[0],
                                    quaternion[1],
                                    quaternion[2],
                                    quaternion[3] );
        positionMsg->Pack();
        socket->Send( positionMsg->GetPackPointer(),
                      positionMsg->GetPackSize() );
        igtl::Sleep(interval); // wait
        }
      }
    }

  //------------------------------------------------------------
  // Close connection 
  if (socket.IsNotNull())
    {
    socket->CloseSocket();
    }

  return EXIT_SUCCESS; 
}

void GenerateRandomOrientation( double orientation[4]  )
{
  // random orientation
  static float theta = 0.0;
  orientation[0]=0.0;
  orientation[1]=0.6666666666*cos(theta);
  orientation[2]=0.577350269189626;
  orientation[3]=0.6666666666*sin(theta);
  theta = theta + 0.1;
}

void GenerateRandomPosition( double position[3] )
{
  //random position
  static float phi = 0.0;
  position[0] = 50.0 * cos(phi);
  position[1] = 50.0 * sin(phi);
  position[2] = 50.0 * cos(phi);
  phi = phi + 2.0;
}
