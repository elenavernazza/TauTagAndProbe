#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <TNtuple.h>
#include <TString.h>
#include <bitset>
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
//#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"

#include "FWCore/Framework/interface/ESHandle.h"  // Handle to read geometry
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"

#include "CondFormats/DataRecord/interface/L1RPCConfigRcd.h"
#include "CondFormats/L1TObjects/interface/L1RPCConfig.h"
#include "CondFormats/DataRecord/interface/L1RPCConeBuilderRcd.h"
#include "CondFormats/RPCObjects/interface/L1RPCConeBuilder.h"
#include "CondFormats/DataRecord/interface/L1RPCHwConfigRcd.h"
#include "CondFormats/RPCObjects/interface/L1RPCHwConfig.h"
#include "CondFormats/DataRecord/interface/L1RPCBxOrConfigRcd.h"
#include "CondFormats/L1TObjects/interface/L1RPCBxOrConfig.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

class PFMetNoMuProducer : public edm::stream::EDProducer<> {
public:
explicit PFMetNoMuProducer(const edm::ParameterSet&);
  ~PFMetNoMuProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  const edm::EDGetTokenT<pat::METCollection> thePFMETCollection_;
  const edm::EDGetTokenT<pat::MuonCollection> theMuonCollection_;
  
private:
  virtual void beginStream(edm::StreamID) override;
virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
PFMetNoMuProducer::PFMetNoMuProducer(const edm::ParameterSet& iConfig) :
thePFMETCollection_ (consumes<pat::METCollection>  (iConfig.getParameter<edm::InputTag>("pfMETCollection"))),
theMuonCollection_  (consumes<pat::MuonCollection> (iConfig.getParameter<edm::InputTag>("muonCollection")))
{
  produces<pat::METCollection>();
  //register your products
  /* Examples
     produces<ExampleData2>();
     //if do put with a label
     produces<ExampleData2>("label");
 
     //if you want to put into the Run
     produces<ExampleData2,InRun>();
  */
  //now do what ever other initialization is needed
  
}


PFMetNoMuProducer::~PFMetNoMuProducer()
{
 
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PFMetNoMuProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<std::vector<pat::MET>> pfMet;
  iEvent.getByToken(thePFMETCollection_, pfMet);
    
  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(theMuonCollection_, muons);
  
  if (!pfMet.isValid()) {
    edm::LogWarning("L1TPFMetNoMuProducer") << "invalid collection for pfMet" << std::endl;
    return;
  }
  if (!muons.isValid()) {
    edm::LogWarning("L1TPFMetNoMuProducer") << "invalid collection for muons" << std::endl;
    return;
  }
    
  pat::MET thePFMetNoMu = pfMet.product()->front();
  double pfMetNoMuPx = thePFMetNoMu.px();
  double pfMetNoMuPy = thePFMetNoMu.py();
   
  double muPx(0.), muPy(0.);
  // std::cout<<"trying muons"<<std::endl;
    
  for (auto muon = muons->begin(); muon != muons->end(); ++muon) {
    if (muon->isPFMuon()) {
      // std::cout<<"there is a muon!"<<std::endl;
      muPx += muon->px();
      muPy += muon->py();
    }
  }
    
  pfMetNoMuPx += muPx;
  pfMetNoMuPy += muPy;
  math::XYZTLorentzVector pfMetNoMuP4(pfMetNoMuPx, pfMetNoMuPy, 0, hypot(pfMetNoMuPx, pfMetNoMuPy));
    
  thePFMetNoMu.setP4(pfMetNoMuP4);
    
  std::unique_ptr<pat::METCollection> product(new pat::METCollection);
  product->emplace_back(thePFMetNoMu);//.getSpecific(), thePFMetNoMu.sumEt(), thePFMetNoMu.p4(), thePFMetNoMu.vertex());
    
  //iEvent.put(std::move(product), "slimmedMETs");
  iEvent.put(std::move(product));
 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
PFMetNoMuProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
PFMetNoMuProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  PFMetNoMuProducer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
  void
  PFMetNoMuProducer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  PFMetNoMuProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  PFMetNoMuProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PFMetNoMuProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFMetNoMuProducer);
