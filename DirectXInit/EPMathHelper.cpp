#include "EPMathHelper.h"
#include <math.h>

using namespace EPMathHelperFunctions;

//copied over from xnamath's defines
const float EPMathHelper::EP_PI				   = 3.141592654f	  ;
const float EPMathHelper::EP_2PI = 6.283185307f					  ;
const float EPMathHelper::EP_1DIVPI = 0.318309886f				  ;
const float EPMathHelper::EP_1DIV2PI = 0.159154943f				  ;
const float EPMathHelper::EP_PIDIV2 = 1.570796327f				  ;
const float EPMathHelper::EP_PIDIV4 = 0.785398163f				  ;

float EPMathHelper::CyclThetaAngleInPlaneFrom(float x, float z)
{

	float arctan = atanf(z / x);
	float retTheta = 0.0f;

	//we can write this a bit differently to eliminate two branches
	if (x < 0.0f)
	{
		//whether it returns a negative angle, meaning quad 2, or a positive angle, meaning quad 3, we still have to add PI
		retTheta = arctan + EP_PI;

	}
	else
	{
		//check for quadrant 4, else quad 1
		if (arctan < 0.0f)
			retTheta = arctan + EP_2PI;
	}

	return retTheta;

}