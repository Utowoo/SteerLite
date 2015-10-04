// Copyright (c) 2015 Mahyar Khayatkhoei
// Copyright (c) 2009-2014 Shawn Singh, Glen Berseth, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//
   
#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"
   
using namespace Util;
   
Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
    controlPoints.push_back(startPoint);
}
   
Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
    controlPoints = inputPoints;
    sortControlPoints();
}
   
// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
    controlPoints.push_back(inputPoint);
    sortControlPoints();
}
   
// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
    for (int i = 0; i < inputPoints.size(); i++)
        controlPoints.push_back(inputPoints[i]);
    sortControlPoints();
}
   
// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
 
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
     
     
#ifdef ENABLE_GUI
     
 
      //td::cout<<"type: "<<type<<std::endl;
    window = 1;
    
    for(int i = 0; i < controlPoints.size() - 1; i ++)
    {
        
    Point outputPoint;
        Point perviousPoint;
        perviousPoint.x = controlPoints[i].position.x;
        perviousPoint.y = controlPoints[i].position.y;
        perviousPoint.z = controlPoints[i].position.z;
               
        int insertNum = (int)(controlPoints[i + 1].time - controlPoints[i].time) / window;
               
        for(int j = 0; j <= insertNum; j++)
        {
            calculatePoint(outputPoint, controlPoints[i].time+j*window);
            DrawLib::drawLine(perviousPoint, outputPoint, curveColor,5);
            //std::cerr <<"(" <<outputPoint.x <<", "<< outputPoint.z <<")"<< std::endl;
            //std::cerr << outputPoint.z << std::endl;
            perviousPoint.x = outputPoint.x;
            perviousPoint.y = outputPoint.y;
            perviousPoint.z = outputPoint.z;
        }
  
    }
    return;
#endif
}
   
// Sort controlPoints vector in ascending order: min-firsts
void Curve::sortControlPoints()
{
    //================DELETE THIS PART AND THEN START CODING===================
    static bool flag = false;
    if (!flag)
    {   
        //std::cerr << controlPoints[0].time << std::endl;        
        //std::cerr << controlPoints[1].time << std::endl;
        //std::cerr << controlPoints[2].time << std::endl;
        //std::cerr << controlPoints[3].time << std::endl;
        //std::cerr << controlPoints.size() << std::endl;
        //std::cerr << "ERROR>>>>Member function sortControlPoints is not implemented!" << std::endl;
        flag = true;
    }
    //=========================================================================
 
    for(int i = 0; i < controlPoints.size(); i ++){
 
        for(int j = i+1; j < controlPoints.size(); j ++){
            Point temp_position;
            float temp_time;
            Vector temp_tangent;
            if(controlPoints[j].time < controlPoints[i].time){
                temp_position = controlPoints[i].position;
                temp_time = controlPoints[i].time;
                temp_tangent = controlPoints[i].tangent;
                controlPoints[i].position = controlPoints[j].position;
                controlPoints[i].tangent = controlPoints[j].tangent;
                controlPoints[i].time = controlPoints[j].time;
                controlPoints[j].position = temp_position;
                controlPoints[j].tangent = temp_tangent;
                controlPoints[j].time = temp_time;
             
             
            }
        }
    }
  
     
    int fl = 0;
    int count [controlPoints.size()];
    int size  = -1;
 
    for(int i = 1; i < controlPoints.size(); i ++){
        if(controlPoints[i].time  == controlPoints[fl].time){
            size ++;
            count[size] = i;
 
        }
        else{
            fl = i;
        }
    }
 
    for(int i = size; i >= 0;i --){
        controlPoints.erase(controlPoints.begin() + count[i]);
    }
    
  
    return;
}
   
// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
bool Curve::calculatePoint(Point& outputPoint, float time)
{
    // Robustness: make sure there is at least two control point: start and end points
    if (!checkRobust())
        return false;
   
    // Define temporary parameters for calculation
    unsigned int nextPoint;
    float normalTime, intervalTime;
   
    // Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
    if (!findTimeInterval(nextPoint, time))
        return false;
   
    // Calculate position at t = time on curve
    if (type == hermiteCurve)
    {
        outputPoint = useHermiteCurve(nextPoint, time);
    }
    else if (type == catmullCurve)
    {
        outputPoint = useCatmullCurve(nextPoint, time);
    }
   
    // Return
    return true;
}
   
// Check Roboustness
bool Curve::checkRobust()
{
   
    if(controlPoints.size() <= 1)
        return false;
   
    return true;
}
   
// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
    for(int i = 0; i < controlPoints.size(); i++){
        if(time >= controlPoints[i].time)
            nextPoint = i+1;    
    }
    //std::cerr << time << "is the time and next point is "<<nextPoint << std::endl;
   
    return true;
}
   
// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
    Point newPosition;
    float normalTime, intervalTime;
    float h1,h2,h3,h4;
    float x1,y1,z1;
    float x2,y2,z2;
    float difx,dify,difz,distance;
      
    int i = 0;
    i = ((int)time/50) * 10;
    x1 = controlPoints[nextPoint-1].position.x;
    y1 = controlPoints[nextPoint-1].position.y;
    z1 = controlPoints[nextPoint-1].position.z;
         
    x2 = controlPoints[nextPoint].position.x;
    y2 = controlPoints[nextPoint].position.y;
    z2 = controlPoints[nextPoint].position.z;
        
    difx = x2 - x1;
    dify = y2 - y1;
    difz = z2 - z1;
        
    distance = 50;
    //sqrt(difx*difx + dify*dify + difz*difz)*50;
          
    intervalTime = ((time - controlPoints[nextPoint-1].time)/(controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time));
    //std::cerr << time << "intervalTime "<<intervalTime << std::endl;
        
    h1 = 2.0 * intervalTime*intervalTime*intervalTime - 3.0*intervalTime*intervalTime + 1.0;
    h2 = intervalTime*intervalTime*intervalTime -2.0*intervalTime*intervalTime + intervalTime;
    h3 = -2.0 * intervalTime*intervalTime*intervalTime + 3.0*intervalTime*intervalTime;
    h4 = intervalTime*intervalTime*intervalTime -intervalTime*intervalTime;
          
    newPosition.x = h1 * x1 + h2 * controlPoints[nextPoint-1].tangent.x*distance  + h3 * x2 + h4 * controlPoints[nextPoint].tangent.x*distance ;
    newPosition.y = h1 * y1 + h2 * controlPoints[nextPoint-1].tangent.y*distance  + h3 * y2 + h4 * controlPoints[nextPoint].tangent.y*distance ;
    newPosition.z = h1 * z1 + h2 * controlPoints[nextPoint-1].tangent.z*distance  + h3 * z2 + h4 * controlPoints[nextPoint].tangent.z*distance ;
       
    return newPosition; 
}
 
 
// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
    //calculate 
    Point newPosition;
    float normalTime, intervalTime;
    float h1,h2,h3,h4;
    float x1,y1,z1;
    float x2,y2,z2;
    float difx,dify,difz,distance;
 
    if(nextPoint >= 1 & nextPoint <= controlPoints.size()-2){
        float x0, x1, x2, y0, y1, y2, z0, z1, z2;
        float t0, t1, t2;
        unsigned int i = nextPoint;
 
        x0 = controlPoints[i-1].position.x; x1 = controlPoints[i].position.x; x2 = controlPoints[i+1].position.x;
        y0 = controlPoints[i-1].position.y; y1 = controlPoints[i].position.y; y2 = controlPoints[i+1].position.y;
        z0 = controlPoints[i-1].position.z; z1 = controlPoints[i].position.z; z2 = controlPoints[i+1].position.z;
        t0 = controlPoints[i-1].time; t1 = controlPoints[i].time; t2 = controlPoints[i+1].time;
 
        if(nextPoint > 1 & nextPoint < controlPoints.size()-2){
            controlPoints[i].tangent.x = (t1 - t0)/(t2 - t0) * (x2 - x1)/(t2 - t1) + (t2 - t1)/(t2 - t0)*(x1 - x0)/(t1 - t0);
            controlPoints[i].tangent.y = (t1 - t0)/(t2 - t0) * (y2 - y1)/(t2 - t1) + (t2 - t1)/(t2 - t0)*(y1 - y0)/(t1 - t0);
            controlPoints[i].tangent.z = (t1 - t0)/(t2 - t0) * (z2 - z1)/(t2 - t1) + (t2 - t1)/(t2 - t0)*(z1 - z0)/(t1 - t0); 
              
        }
        else if(nextPoint == 1){   
 
            controlPoints[i].tangent.x = (t1 - t0)/(t2 - t0) * (x2 - x1)/(t2 - t1) + (t2 - t1)/(t2 - t0)*(x1 - x0)/(t1 - t0);
            controlPoints[i].tangent.y = (t1 - t0)/(t2 - t0) * (y2 - y1)/(t2 - t1) + (t2 - t1)/(t2 - t0)*(y1 - y0)/(t1 - t0);
            controlPoints[i].tangent.z = (t1 - t0)/(t2 - t0) * (z2 - z1)/(t2 - t1) + (t2 - t1)/(t2 - t0)*(z1 - z0)/(t1 - t0); 
 
            controlPoints[i-1].tangent.x = (t2 - t0)/(t2 - t1) * (x1 - x0)/(t1 - t0) - (t1 - t0)/(t2 - t1)*(x2 - x0)/(t2 - t0);
            controlPoints[i-1].tangent.y = (t2 - t0)/(t2 - t1) * (y1 - y0)/(t1 - t0) - (t1 - t0)/(t2 - t1)*(y2 - y0)/(t2 - t0);
            controlPoints[i-1].tangent.z = (t2 - t0)/(t2 - t1) * (z1 - z0)/(t1 - t0) - (t1 - t0)/(t2 - t1)*(z2 - z0)/(t2 - t0);
  
        }
        else{
 
            controlPoints[i].tangent.x = (t1 - t0)/(t2 - t0) * (x2 - x1)/(t2 - t1) + (t2 - t1)/(t2 - t0)*(x1 - x0)/(t1 - t0);
            controlPoints[i].tangent.y = (t1 - t0)/(t2 - t0) * (y2 - y1)/(t2 - t1) + (t2 - t1)/(t2 - t0)*(y1 - y0)/(t1 - t0);
            controlPoints[i].tangent.z = (t1 - t0)/(t2 - t0) * (z2 - z1)/(t2 - t1) + (t2 - t1)/(t2 - t0)*(z1 - z0)/(t1 - t0); 
            
            controlPoints[i+1].tangent.x = (t0 - t2)/(t0 - t1) * (x1 - x2)/(t1 - t2) - (t1 - t2)/(t0 - t1)*(x0 - x2)/(t0 - t2);
            controlPoints[i+1].tangent.y = (t0 - t2)/(t0 - t1) * (y1 - y2)/(t1 - t2) - (t1 - t2)/(t0 - t1)*(y0 - y2)/(t0 - t2);
            controlPoints[i+1].tangent.z = (t0 - t2)/(t0 - t1) * (z1 - z2)/(t1 - t2) - (t1 - t2)/(t0 - t1)*(z0 - z2)/(t0 - t2);
             
        }
 
  }
     
  
     
    //previous point
    x1 = controlPoints[nextPoint-1].position.x;
    y1 = controlPoints[nextPoint-1].position.y;
    z1 = controlPoints[nextPoint-1].position.z;
     
    //current point
    x2 = controlPoints[nextPoint].position.x;
    y2 = controlPoints[nextPoint].position.y;
    z2 = controlPoints[nextPoint].position.z;
     
    // to calculate distance
    difx = x2 - x1;
    dify = y2 - y1;
    difz = z2 - z1;
        
    distance = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
    // sqrt(difx*difx + dify*dify + difz*difz) ;
 
    // to calculate t 
    intervalTime = ((time - controlPoints[nextPoint-1].time)/(controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time));
    //std::cerr << time << "intervalTime "<<intervalTime << std::endl;
 
    // to calculate tangent of previous point and current point
     
 
    h1 = 2.0 * intervalTime*intervalTime*intervalTime - 3.0*intervalTime*intervalTime + 1.0;
    h2 = intervalTime*intervalTime*intervalTime -2.0*intervalTime*intervalTime + intervalTime;
    h3 = -2.0 * intervalTime*intervalTime*intervalTime + 3.0*intervalTime*intervalTime;
    h4 = intervalTime*intervalTime*intervalTime -intervalTime*intervalTime;
          
    newPosition.x = h1 * x1 + h2 * controlPoints[nextPoint-1].tangent.x*distance  + h3 * x2 + h4 * controlPoints[nextPoint].tangent.x*distance ;
    newPosition.y = h1 * y1 + h2 * controlPoints[nextPoint-1].tangent.y*distance  + h3 * y2 + h4 * controlPoints[nextPoint].tangent.y*distance ;
    newPosition.z = h1 * z1 + h2 * controlPoints[nextPoint-1].tangent.z*distance  + h3 * z2 + h4 * controlPoints[nextPoint].tangent.z*distance ;
    
    return newPosition; 
 
 
}
