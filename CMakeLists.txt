CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

PROJECT(vtkImageNormalizedCrossCorrelation)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(vtkImageNormalizedCrossCorrelationExample vtkImageNormalizedCrossCorrelationExample.cxx
vtkImageNormalizedCrossCorrelation.cxx)
TARGET_LINK_LIBRARIES(vtkImageNormalizedCrossCorrelationExample vtkHybrid)

