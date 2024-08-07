#include "NNClusteringProcessor.hh"
#include <iostream>

#include "time.h"

#include <IMPL/LCCollectionVec.h>
#include <IMPL/ClusterImpl.h>

#include "NNClusters.h"

#include <algorithm>

using namespace lcio ;
using namespace marlin ;

NNClusteringProcessor aNNClusteringProcessor ;


NNClusteringProcessor::NNClusteringProcessor() : Processor("NNClusteringProcessor") {
  
  // modify processor description
  _description = "NNClusteringProcessor : simple nearest neighbour clustering" ;
  

  StringVec colDefault ;
  colDefault.push_back("ecal" ) ;
  
  registerInputCollections( LCIO::CALORIMETERHIT,
			    "EcalCollection" , 
			    "Name of the input collections"  ,
			    _colNames ,
			    colDefault ) ;
  
  registerOutputCollection( LCIO::CLUSTER,
			    "OutputCollection" , 
			    "Name of the output collections"  ,
			    _outputColName ,
			    std::string("NNClusters" ) ) ;
 

  registerProcessorParameter( "DistanceCut" , 
			      "Cut for distance between hits in mm"  ,
			      _distCut ,
			       (float) 40.0 ) ;

  registerProcessorParameter( "EnergyCut" , 
			      "Cut for hit energy in GeV"  ,
			      _eCut ,
			       (float) 0.0 ) ;

}


void NNClusteringProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void NNClusteringProcessor::processRunHeader( LCRunHeader*  /*run*/) { 

  _nRun++ ;
} 

void NNClusteringProcessor::processEvent( LCEvent * evt ) { 


  std::cout << " ---- NNClusteringProcessor::processEvent() - evt: " 
	    << evt->getRunNumber() << " , " << evt->getEventNumber() 
	    << std::endl ;

  clock_t start =  clock () ; 


  LCCollectionVec* lcioClusters = new LCCollectionVec( LCIO::CLUSTER )  ;
  
  GenericHitVec<CalorimeterHit> h ;

  GenericClusterVec<CalorimeterHit> cl ;
  
  EnergyCut<CalorimeterHit> eCut( _eCut ) ;
  
  ZIndex<CalorimeterHit,100> zIndex( 0. , 250. ) ; // this is modified wrt the Marlin Reco -4300,4300

  NNDistance< CalorimeterHit, float> dist( _distCut )  ;

  LCIOCluster<CalorimeterHit> converter ;
  

  int nHit = 0 ;
  // create a vector of generic hits from the collection applying an energy cut
  for( StringVec::iterator it = _colNames.begin() ; it !=  _colNames.end() ; it++ ){  
    
    LCCollection* col =  evt->getCollection( *it )  ;
    nHit += col->getNumberOfElements()  ;

//     addToGenericHitVec( h , col , eCut ) ;
    addToGenericHitVec( h , col , eCut , zIndex ) ;
  }



  // cluster the hits with a nearest neighbour condition
  cluster( h.begin() , h.end() , std::back_inserter( cl )  , &dist ) ;
  
  std::cout << "  passing " << h.size() << " of " << nHit  
	    << "  hits to clustering (E_cut: " << _eCut << ") " 
	    << "  found  " << cl.size() << " clusters " << std::endl ;

  // create lcio::Clusters from the clustered GenericHits
  std::transform( cl.begin(), cl.end(), std::back_inserter( *lcioClusters ) , converter ) ;


  evt->addCollection( lcioClusters , _outputColName ) ;
  
  
  _nEvt ++ ;

  clock_t end =  clock () ; 
  
  std::cout << " ---- NNClusterProcessor::processEvent() - time: " 
	    <<  double( end - start ) / double(CLOCKS_PER_SEC) 
	    << std::endl  ;

}



void NNClusteringProcessor::check( LCEvent *  /*evt*/ ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void NNClusteringProcessor::end(){ 
  
//   std::cout << "NNClusterProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}