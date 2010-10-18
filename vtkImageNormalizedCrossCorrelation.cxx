/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageNormalizedCrossCorrelation.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkImageNormalizedCrossCorrelation.h"

#include "vtkSmartPointer.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkImageCast.h"
#include "vtkImageExtractComponents.h"
#include "vtkImageAppendComponents.h"
#include "vtkImageShiftScale.h"
#include "vtkImageMathematics.h"

#include <vector>

vtkStandardNewMacro(vtkImageNormalizedCrossCorrelation);

vtkImageNormalizedCrossCorrelation::vtkImageNormalizedCrossCorrelation()
{
  this->SetNumberOfInputPorts(2);
}

int vtkImageNormalizedCrossCorrelation::RequestData(vtkInformation *vtkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector)
{
  // Get the info objects
  vtkInformation *image1Info = inputVector[0]->GetInformationObject(0);
  vtkInformation *image2Info = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Get the input and ouptut
  vtkImageData *input1 = vtkImageData::SafeDownCast(
      image1Info->Get(vtkDataObject::DATA_OBJECT()));

  vtkImageData *input2 = vtkImageData::SafeDownCast(
      image2Info->Get(vtkDataObject::DATA_OBJECT()));

  vtkImageData *output = vtkImageData::SafeDownCast(
      outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Cast both images to double
  vtkSmartPointer<vtkImageCast> image1CastFilter =
    vtkSmartPointer<vtkImageCast>::New();
  image1CastFilter->SetInputConnection(input1->GetProducerPort());
  image1CastFilter->SetOutputScalarTypeToDouble();
  image1CastFilter->Update();

  vtkSmartPointer<vtkImageData> image1 =
    vtkSmartPointer<vtkImageData>::New();
  image1->DeepCopy(image1CastFilter->GetOutput());

  vtkSmartPointer<vtkImageCast> image2CastFilter =
    vtkSmartPointer<vtkImageCast>::New();
  image2CastFilter->SetInputConnection(input2->GetProducerPort());
  image2CastFilter->SetOutputScalarTypeToDouble();
  image2CastFilter->Update();

  vtkSmartPointer<vtkImageData> image2 =
    vtkSmartPointer<vtkImageData>::New();
  image2->DeepCopy(image2CastFilter->GetOutput());

  // http://www.mathworks.com/help/toolbox/images/ref/normxcorr2.html

  // Compute and subtract the mean from both images
  SubtractMean(image1);
  SubtractMean(image2);

  return 1;
}

double vtkImageNormalizedCrossCorrelation::SubtractMean(vtkImageData* image)
{
  std::vector<vtkSmartPointer<vtkImageData> > normalizedComponents;

  for(vtkIdType c = 0; c < image->GetNumberOfScalarComponents(); c++)
    {
    vtkSmartPointer<vtkImageExtractComponents> extractComponentFilter =
      vtkSmartPointer<vtkImageExtractComponents>::New();
    extractComponentFilter->SetInputConnection(image->GetProducerPort());
    extractComponentFilter->SetComponents(c);
    extractComponentFilter->Update();

    vtkSmartPointer<vtkImageData> component =
      vtkSmartPointer<vtkImageData>::New();
    component->DeepCopy(extractComponentFilter->GetOutput());

    double mean = ComputeMean(component);

    vtkSmartPointer<vtkImageShiftScale> shiftScaleFilter =
      vtkSmartPointer<vtkImageShiftScale>::New();
    shiftScaleFilter->SetOutputScalarTypeToUnsignedChar();
    shiftScaleFilter->SetInputConnection(component->GetProducerPort());
    shiftScaleFilter->SetShift(-mean);
    shiftScaleFilter->Update();

    vtkSmartPointer<vtkImageData> shiftedComponent =
      vtkSmartPointer<vtkImageData>::New();
    shiftedComponent->DeepCopy(shiftScaleFilter->GetOutput());

    normalizedComponents.push_back(shiftedComponent);
    }

  vtkSmartPointer<vtkImageAppendComponents> appendFilter =
    vtkSmartPointer<vtkImageAppendComponents>::New();
  for(vtkIdType c = 0; c < image->GetNumberOfScalarComponents(); c++)
    {
    appendFilter->AddInputConnection(0, normalizedComponents[c]->GetProducerPort());
    }
  appendFilter->Update();

  image->DeepCopy(appendFilter->GetOutput());
}

double vtkImageNormalizedCrossCorrelation::ComputeMean(vtkImageData* image)
{
  // This function expects the image type to be double and to have 1 scalar component
  int extent[6];
  image->GetExtent(extent);

  double pixelSum = 0;
  unsigned int pixelCount = 0;

  for(unsigned int i = extent[0]; i <= extent[1]; i++)
    {
    for(unsigned int j = extent[2]; j <= extent[3]; j++)
      {
      for(unsigned int k = extent[4]; k <= extent[5]; k++)
        {
        double* pixel = static_cast<double*>(image->GetScalarPointer(i,j,k));
        pixelSum += pixel[0];
        pixelCount++;
        }
      }
    }

  return pixelSum/static_cast<double>(pixelCount);
}

void vtkImageNormalizedCrossCorrelation::CrossCorrelation(vtkImageData* image, vtkImageData* patch, vtkImageData* output)
{
  // http://www.mathworks.com/help/toolbox/images/ref/normxcorr2.html

  vtkSmartPointer<vtkImageData> normalizedImage =
    vtkSmartPointer<vtkImageData>::New();
  NormalizeImage(image, normalizedImage);

  vtkSmartPointer<vtkImageData> normalizedPatch =
    vtkSmartPointer<vtkImageData>::New();
  NormalizeImage(patch, normalizedPatch);

  vtkSmartPointer<vtkImageMathematics> squareImage =
    vtkSmartPointer<vtkImageMathematics>::New();
  squareImage->SetOperationToSquare();
  squareImage->SetInput(normalizedImage);
  squareImage->Update();

  vtkSmartPointer<vtkImageMathematics> squarePatch =
    vtkSmartPointer<vtkImageMathematics>::New();
  squarePatch->SetOperationToSquare();
  squarePatch->SetInput(normalizedPatch);
  squarePatch->Update();
}

void vtkImageNormalizedCrossCorrelation::NormalizeImage(vtkImageData* input, vtkImageData* output)
{
  // Cast both images to double
  vtkSmartPointer<vtkImageCast> imageCastFilter =
    vtkSmartPointer<vtkImageCast>::New();
  imageCastFilter->SetInputConnection(input->GetProducerPort());
  imageCastFilter->SetOutputScalarTypeToDouble();
  imageCastFilter->Update();

  vtkSmartPointer<vtkImageData> doubleImage =
    vtkSmartPointer<vtkImageData>::New();
  doubleImage->DeepCopy(imageCastFilter->GetOutput());

  // Compute and subtract the mean
  SubtractMean(doubleImage);

  output->ShallowCopy(doubleImage);
}

double vtkImageNormalizedCrossCorrelation::PixelSum(vtkImageData* image)
{
  double pixelSum = 0;

  int extent[6];
  image->GetExtent(extent);

  for(unsigned int i = extent[0]; i <= extent[1]; i++)
    {
    for(unsigned int j = extent[2]; j <= extent[3]; j++)
      {
      for(unsigned int k = extent[4]; k <= extent[5]; k++)
        {
        double* pixel = static_cast<double*>(image->GetScalarPointer(i,j,k));
        pixelSum += pixel[0];
        }
      }
    }

  return pixelSum;
}

void vtkImageNormalizedCrossCorrelation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}