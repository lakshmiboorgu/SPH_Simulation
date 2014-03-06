
#include "SPHNode.h"

MTypeId    SPHNode::id( 0x00000231 );
MObject     SPHNode::aOutValue;
MObject    SPHNode::aInValue;
MObject     SPHNode::aMagnitude;
MObject     SPHNode::aMean;
MObject     SPHNode::aVariance;

SPHNode::SPHNode() 
{
}


SPHNode::~SPHNode() 
{
}


void* SPHNode::creator()
{ 
    return new SPHNode(); 
}


MStatus SPHNode::compute( const MPlug& plug, MDataBlock& data )
{
    MStatus status;

    if ( plug != aOutValue )
    {
        return MS::kUnknownParameter;
    }


    float inputValue = data.inputValue( aInValue, &status ).asFloat();
    float magnitude = data.inputValue( aMagnitude, &status ).asFloat();
    float mean = data.inputValue( aMean, &status ).asFloat();
    float variance = data.inputValue( aVariance, &status ).asFloat();
    if ( variance <= 0.0f )
    {
        variance = 0.001f;
    }

    float xMinusB = inputValue - mean;
    float output = exp( -(xMinusB * xMinusB) / (2.0f * variance) );

	

	//std::cout<<"output valuse is :"<< output;

    MDataHandle hOutput = data.outputValue( aOutValue, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    hOutput.setFloat( output );
    hOutput.setClean();
    data.setClean( plug );

    return MS::kSuccess;
}


MStatus SPHNode::initialize()
{
    MStatus status;
    MFnNumericAttribute nAttr;

    aOutValue = nAttr.create( "outValue", "outValue", MFnNumericData::kFloat );
    nAttr.setWritable( false );
    nAttr.setStorable( false );
    addAttribute( aOutValue );

    aInValue = nAttr.create( "inValue", "inValue", MFnNumericData::kFloat );
    nAttr.setKeyable( true );
    addAttribute( aInValue );
    attributeAffects( aInValue, aOutValue );

    aMagnitude = nAttr.create( "magnitude", "magnitude", MFnNumericData::kFloat );
    nAttr.setKeyable( true );
    addAttribute( aMagnitude );
    attributeAffects( aMagnitude, aOutValue );

    aMean = nAttr.create( "mean", "mean", MFnNumericData::kFloat );
    nAttr.setKeyable( true );
    addAttribute( aMean );
    attributeAffects( aMean, aOutValue );

    aVariance = nAttr.create( "variance", "variance", MFnNumericData::kFloat );
    nAttr.setKeyable( true );
    addAttribute( aVariance );
    attributeAffects( aVariance, aOutValue );

    return MS::kSuccess;
}
