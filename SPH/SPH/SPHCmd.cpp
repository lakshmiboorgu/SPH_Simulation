#include "SPHCmd.h"


SPHCommand::SPHCommand()
{
}


void* SPHCommand::creator()
{
    return new SPHCommand;
}


MStatus SPHCommand::doIt( const MArgList& argList )
{
    MGlobal::displayInfo( "I just wrote my first command plug-in!" );
    return MS::kSuccess;
}