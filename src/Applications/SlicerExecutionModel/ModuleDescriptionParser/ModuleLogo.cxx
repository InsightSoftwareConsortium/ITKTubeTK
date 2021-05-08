/*=========================================================================

  Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   Module Description Parser
  Module:    $HeadURL: http://svn.slicer.org/Slicer4/trunk/Libs/ModuleDescriptionParser/ModuleDescription.h $
  Date:      $Date: 2006-11-10 15:42:21 -0500 (Fri, 10 Nov 2006) $
  Version:   $Revision: 1562 $

==========================================================================*/

#include "ModuleLogo.h"

//----------------------------------------------------------------------------
ModuleLogo::ModuleLogo()
  : Width(0), Height(0), PixelSize(0), BufferLength(0), Options(0), Logo("")
{
}

//----------------------------------------------------------------------------
ModuleLogo::ModuleLogo(const ModuleLogo& logo)
{
  this->Width = logo.Width;
  this->Height = logo.Height;
  this->PixelSize = logo.PixelSize;
  this->BufferLength = logo.BufferLength;
  this->Options = logo.Options;
  this->Logo = logo.Logo;
}

//----------------------------------------------------------------------------
void
ModuleLogo::operator=(const ModuleLogo &logo)
{
  this->Width = logo.Width;
  this->Height = logo.Height;
  this->PixelSize = logo.PixelSize;
  this->BufferLength = logo.BufferLength;
  this->Options = logo.Options;
  this->Logo = logo.Logo;
}

//----------------------------------------------------------------------------
ModuleLogo::~ModuleLogo()
{

}

//----------------------------------------------------------------------------
void
ModuleLogo
::SetLogo(char const * logo, int width, int height, int pixelSize, unsigned long bufferLength, int options)
{
  this->Width = width;
  this->Height = height;
  this->PixelSize = pixelSize;
  this->BufferLength = bufferLength;
  this->Options = options;
  this->Logo = std::string(logo, bufferLength);
}

//----------------------------------------------------------------------------
int
ModuleLogo
::GetWidth() const
{
  return this->Width;
}

//----------------------------------------------------------------------------
int
ModuleLogo
::GetHeight() const
{
  return this->Height;
}

//----------------------------------------------------------------------------
int
ModuleLogo
::GetPixelSize() const
{
  return this->PixelSize;
}

//----------------------------------------------------------------------------
unsigned long
ModuleLogo
::GetBufferLength() const
{
  return this->BufferLength;
}

//----------------------------------------------------------------------------
int
ModuleLogo
::GetOptions() const
{
  return this->Options;
}

//----------------------------------------------------------------------------
const char *
ModuleLogo
::GetLogo() const
{
  return this->Logo.c_str();
}
