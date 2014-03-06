

#include <maya/MPxCommand.h>
#include <maya/MGlobal.h>
#include <maya/MObject.h>

class SPHCommand : public MPxCommand
{
public:
    SPHCommand();
    virtual MStatus doIt( const MArgList& argList );
    static void* creator();
};


