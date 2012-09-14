// -*- C++ -*-
//
// Package:    RutgersSubJetValidator
// Class:      RutgersSubJetValidator
// 
/**\class RutgersSubJetValidator RutgersSubJetValidator.cc RutgersSandbox/RutgersSubJetValidator/src/RutgersSubJetValidator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dinko Ferencek
//         Created:  Mon Aug 13 11:13:25 CDT 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

// ROOT
#include <TLorentzVector.h>

//
// class declaration
//

class RutgersSubJetValidator : public edm::EDAnalyzer {
   public:
      explicit RutgersSubJetValidator(const edm::ParameterSet&);
      ~RutgersSubJetValidator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

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
RutgersSubJetValidator::RutgersSubJetValidator(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


RutgersSubJetValidator::~RutgersSubJetValidator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RutgersSubJetValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::BasicJetCollection> genJets;
  iEvent.getByLabel(edm::InputTag("ak5GenJetsRU"), genJets);

  edm::Handle<reco::GenJetCollection> genSubJets;
  iEvent.getByLabel(edm::InputTag("ak5GenJetsRU:SubJets"), genSubJets);

  unsigned jet_count = 0;
  reco::BasicJetCollection::const_iterator genJetIt = genJets->begin();
  for( ; genJetIt != genJets->end(); ++genJetIt )
  {
    // exit from loop when you reach the required number of jets
    if(jet_count > 1)
      break;

    unsigned const_size = genJetIt->getJetConstituents().size();
    std::cout << "Jet " << jet_count << ": eta = " << genJetIt->eta() << " phi = " << genJetIt->phi() << " Pt = " << genJetIt->pt() << " Area = " << genJetIt->jetArea() << std::endl
              << "       Constituents = " << const_size << std::endl;

//     TLorentzVector vSum;
//     for(unsigned i = 0; i<const_size; ++i)
//     {
//       std::cout << "  Constituent " << i << ": eta = " << genJetIt->getJetConstituents()[i]->eta() << " phi = " << genJetIt->getJetConstituents()[i]->phi()
//                 << " Pt = " << genJetIt->getJetConstituents()[i]->pt() << std::endl;
//       TLorentzVector v;
//       v.SetPtEtaPhiE(genJetIt->getJetConstituents()[i]->pt(), genJetIt->getJetConstituents()[i]->eta(),
// 		     genJetIt->getJetConstituents()[i]->phi(), genJetIt->getJetConstituents()[i]->energy());
//       vSum += v;
//     }
//     std::cout << "  Constituent sum: eta = " << vSum.Eta() << " phi = " << vSum.Phi() << " Pt = " << vSum.Pt() << std::endl;

    ++jet_count;
  }

  unsigned subjet_count = 0;
  reco::GenJetCollection::const_iterator genSubJetIt = genSubJets->begin();
  for( ; genSubJetIt != genSubJets->end(); ++genSubJetIt )
  {
    // exit from loop when you reach the required number of jets
    if(subjet_count > 3)
      break;

    unsigned const_size = genSubJetIt->getJetConstituents().size();
    std::cout << "SubJet " << subjet_count << ": eta = " << genSubJetIt->eta() << " phi = " << genSubJetIt->phi() << " Pt = " << genSubJetIt->pt() << " Area = " << genSubJetIt->jetArea() << std::endl
              << "       Constituents = " << const_size << std::endl;

    ++subjet_count;
  }

}

// ------------ method called once each job just before starting event loop  ------------
void 
RutgersSubJetValidator::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RutgersSubJetValidator::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
RutgersSubJetValidator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
RutgersSubJetValidator::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
RutgersSubJetValidator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
RutgersSubJetValidator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RutgersSubJetValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RutgersSubJetValidator);
