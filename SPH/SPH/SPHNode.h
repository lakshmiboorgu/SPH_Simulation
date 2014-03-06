
#include <maya/MPxNode.h>
#include <maya/MFnNumericAttribute.h>

#include <math.h>
#include <iostream>


class SPHNode : public MPxNode
{
public:
						SPHNode();
	virtual				~SPHNode(); 
	static  void*		creator();

	virtual MStatus     compute( const MPlug& plug, MDataBlock& data );
	static  MStatus		initialize();

	static MTypeId	id;
    static MObject  aOutValue;
    static MObject  aInValue;
    static MObject  aMagnitude;
    static MObject  aMean;
    static MObject  aVariance;

};


