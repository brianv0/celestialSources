
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"

#include "FluxSvc/IRegisterSource.h"
#include "FluxSvc/ISpectrumFactory.h"
#include "FluxSvc/IFluxSvc.h"
// GRB includes

#include "GRB/GRBmanager.h"
#include "GRBmaker/GRBobsSpectrum.h"


/** @class RegisterGRB
 *  @brief Register the various GRB sources
 *  
 *   @author Toby Burnett
 *   $Header$
 */
class RegisterGRB : public AlgTool, virtual public IRegisterSource {
public:
    
    RegisterGRB( const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~RegisterGRB() { }
    
    /// implement to define sources: will be called from FluxSvc
    StatusCode registerMe(IFluxSvc* );
    
};


// Static factory for instantiation of algtool objects
static ToolFactory<RegisterGRB> s_factory;
const IToolFactory& RegisterGRBFactory = s_factory;

// Standard Constructor
RegisterGRB::RegisterGRB(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : AlgTool( type, name, parent ) {
    
    // Declare additional interface
    declareInterface<IRegisterSource>(this);
        
}


StatusCode RegisterGRB::registerMe(IFluxSvc* fsvc) 
{
    MsgStream  log(msgSvc(), name());

    log << MSG::INFO << "Register GRB Manager..." << endreq;
    static RemoteSpectrumFactory<GRBmanager> factory1(fsvc);
    
    log << MSG::INFO << "Register GRB Obs Spectra..." << endreq;
    static RemoteSpectrumFactory<GRBobsSpectrum> factory2(fsvc);
    
    return StatusCode::SUCCESS;
} 
