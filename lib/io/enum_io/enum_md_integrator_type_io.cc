#include "enum_md_integrator_type_io.h"

#include <string>
using namespace std;
using namespace Chroma;

namespace Chroma { 
  
  namespace MDIntegratorTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theMDIntegratorTypeMap::Instance().registerPair(string("PQP_LEAPFROG"), MD_PQP_LEAPFROG);
      success &=theMDIntegratorTypeMap::Instance().registerPair(string("QPQ_LEAPFROG"), MD_QPQ_LEAPFROG);
      
      return success;
    }

    const string typeIDString = "MDIntegratorType";
    bool registered = registerAll();
  };

  using namespace MDIntegratorTypeEnv;
  //! Read an MDIntegratorType enum
  void read(XMLReader& xml_in,  const string& path, MDIntegratorType& t) {
    theMDIntegratorTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an MDIntegratorType enum
  void write(XMLWriter& xml_out, const string& path, const MDIntegratorType& t) {
    theMDIntegratorTypeMap::Instance().write(typeIDString,xml_out, path, t);
  }
};
