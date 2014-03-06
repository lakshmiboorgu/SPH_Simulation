#include "SPHCmd.h"
#include "SPHNode.h"


#include <maya/MFnPlugin.h>

MStatus initializePlugin( MObject obj )
{
    MStatus status;

    MFnPlugin fnPlugin( obj, "lakshmi boorgu", "1.0", "Any" );

    status = fnPlugin.registerCommand( "myFancyCommand", SPHCommand::creator );
    CHECK_MSTATUS_AND_RETURN_IT( status );

	 status = fnPlugin.registerNode( "SPH",
        SPHNode::id,
        SPHNode::creator,
        SPHNode::initialize );
    CHECK_MSTATUS_AND_RETURN_IT( status );


    return MS::kSuccess;
}


MStatus uninitializePlugin( MObject obj )
{
    MStatus status;

    MFnPlugin fnPlugin( obj );

    status = fnPlugin.deregisterCommand( "myFancyCommand" );
    CHECK_MSTATUS_AND_RETURN_IT( status );

	status = fnPlugin.deregisterNode( SPHNode::id );
    CHECK_MSTATUS_AND_RETURN_IT( status );

    return MS::kSuccess;
}