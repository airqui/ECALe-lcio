#include "NNClusteringProcessor.hh"
#include <iostream>
#include "WeightedPoints3D.h"

#include "time.h"

#include <IMPL/LCCollectionVec.h>
#include <IMPL/ClusterImpl.h>

#include "NNClustersAIQ.h"

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
  
  ZIndex<CalorimeterHit,100> zIndex( 0. , 10. ) ; // this is modified wrt the Marlin Reco -4300,4300

  NNDistance< CalorimeterHit, float> dist( _distCut , _eCut)  ;

  LCIOCluster<CalorimeterHit> converter ;
  

  int nHit = 0 ;
  // create a vector of generic hits from the collection applying an energy cut
  for( StringVec::iterator it = _colNames.begin() ; it !=  _colNames.end() ; it++ ){  
    
    LCCollection* col =  evt->getCollection( *it )  ;
    nHit += col->getNumberOfElements()  ;

//     addToGenericHitVec( h , col , eCut ) ;
    addToGenericHitVec( h , col , eCut);//, zIndex ) ;
  }


  int npoints_index=0;
  int sum_wgt_index=0 ;
  int sum_wgt2_index=0 ;
  int sum_wgt4_index=0 ;
  StringVec shape_keys ;
  // cluster the hits with a nearest neighbour condition
  cluster( h.begin() , h.end() , std::back_inserter( cl )  , &dist ) ;
  
  std::cout << "  passing " << h.size() << " of " << nHit  
	    << "  hits to clustering (E_cut: " << _eCut << ") " 
	    << "  found  " << cl.size() << " clusters " << std::endl ;

  // create lcio::Clusters from the clustered GenericHits
  std::transform( cl.begin(), cl.end(), std::back_inserter( *lcioClusters ) , converter ) ;

  int n_el_col=lcioClusters->getNumberOfElements();

  StringVec shapeParams ;
  shapeParams = lcioClusters->getParameters().getStringVals ("ClusterShapeParameters", shapeParams);
  shapeParams.push_back("npoints") ;
  shapeParams.push_back("sum_wgt") ;
  shapeParams.push_back("sum_wgt^2") ;
  shapeParams.push_back("sum_wgt^4") ;
  lcioClusters->parameters().setValues( "ClusterShapeParameters" , shapeParams ) ;
  shape_keys = lcioClusters->getParameters().getStringVals( std::string("ClusterShapeParameters"),shape_keys);
  for ( unsigned kkk=0 ; kkk < shape_keys.size() ; kkk++ ) {
	  if ( shape_keys[kkk] == "npoints" )   { npoints_index  = kkk ; }
	  if ( shape_keys[kkk] == "sum_wgt" ) { sum_wgt_index = kkk ; }
	  if ( shape_keys[kkk] == "sum_wgt^2" ) { sum_wgt2_index = kkk ; }
    if ( shape_keys[kkk] == "sum_wgt^4" ) { sum_wgt4_index = kkk ; }
  }

  for(int iel=0; iel<n_el_col; iel++) {
    ClusterImpl* cluster_ = dynamic_cast<ClusterImpl*>(lcioClusters->getElementAt(iel));

	  unsigned int n = cluster_->getCalorimeterHits().size();
	    
    std::vector<double> ehit,xhit,yhit,zhit;
    for (int jj=0; jj< n ; jj++){
      CalorimeterHit* calo_hit  = dynamic_cast<CalorimeterHit*>(cluster_->getCalorimeterHits().at(jj));
      ehit.push_back(calo_hit->getEnergy());
      xhit.push_back(calo_hit->getPosition()[0]);
      yhit.push_back(calo_hit->getPosition()[1]);
      zhit.push_back(calo_hit->getPosition()[2]);
    }
    
    WeightedPoints3D wgtp = WeightedPoints3D( n , &ehit[0] , &xhit[0], &yhit[0], &zhit[0] ); 
    double* cog=wgtp.getCentreOfGravity();

    double *covv=wgtp.getCentreOfGravityErrors();
    double cov[3][3]; for( int iii=0 ; iii<3 ; iii++ ) { for ( int jjj=0 ; jjj<3 ; jjj++ ) { cov[iii][jjj]=covv[jjj+iii*3] ; } } ; 

    double *eval=wgtp.getEigenVal(); 
    double *eval_err=wgtp.getEigenValErrors(); 

    double *evpv=wgtp.getEigenVecPolar(); 
	  double evp[2][3]; for( int iii=0 ; iii<2 ; iii++ ) { for ( int jjj=0 ; jjj<3 ; jjj++ ) { evp[iii][jjj]=evpv[jjj+iii*3] ; } } ;

    double *evpc=wgtp.getEigenVecCartesian(); 
    double evc[2][3]; for( int iii=0 ; iii<2 ; iii++ ) { for ( int jjj=0 ; jjj<3 ; jjj++ ) { evc[iii][jjj]=evpc[jjj+iii*3] ; } } ;

    double* evpev=wgtp.getEigenVecPolarErrors();              
    double evpe[2][2][3]; for( int iii=0 ; iii<2 ; iii++ ) { for ( int jjj=0 ; jjj<2 ; jjj++ ) {for ( int kkk=0 ; kkk<3 ; kkk++ ) { evpe[iii][jjj][kkk]=evpev[kkk+jjj*3+iii*3*2] ; } } };

    double sum_wgtsqr = wgtp.getTotalSquaredWeight();
    double sum_wgt4 = wgtp.getTotalQuarticWeight();
    double sum_wgt = wgtp.getTotalWeight();
    
    FloatVec PositionError ;
    PositionError.push_back(cov[0][0]*sum_wgtsqr/(sum_wgt*sum_wgt));
  	PositionError.push_back(cov[0][1]*sum_wgtsqr/(sum_wgt*sum_wgt));
	  PositionError.push_back(cov[1][1]*sum_wgtsqr/(sum_wgt*sum_wgt));
	  PositionError.push_back(cov[0][2]*sum_wgtsqr/(sum_wgt*sum_wgt));
	  PositionError.push_back(cov[1][2]*sum_wgtsqr/(sum_wgt*sum_wgt));
	  PositionError.push_back(cov[2][2]*sum_wgtsqr/(sum_wgt*sum_wgt));
    streamlog_out(DEBUG3) << "i PositionError : "  << std::endl;  streamlog_out(DEBUG3) << "i " ;
    for ( int iii=0 ; iii < 6 ; iii++ ) { streamlog_out(DEBUG3) << PositionError[iii] << " " ;}  ;streamlog_out(DEBUG3) << std::endl;

    float Position[3];
    Position[0] = cog[0];
    Position[1] = cog[1];
    Position[2] = cog[2];
    streamlog_out(DEBUG3) << "i Position : "  << std::endl;  streamlog_out(DEBUG3) << "i " ;
    for ( int iii=0 ; iii < 3 ; iii++ ) { streamlog_out(DEBUG3) << Position[iii] << " " ;}  ;streamlog_out(DEBUG3) << std::endl;

    float theta = evp[0][2];
    float phi = evp[1][2];
    streamlog_out(DEBUG3) << "i theta/phi : "  << std::endl;  streamlog_out(DEBUG3) << "i " ;
    streamlog_out(DEBUG3) << theta << " "  << phi <<  std::endl;

    FloatVec DirectionError ;
	  DirectionError.push_back(evpe[0][0][2]);
	  DirectionError.push_back(evpe[1][0][2]);
	  DirectionError.push_back(evpe[1][1][2]);
    streamlog_out(DEBUG3) << "i DirectionError : "  << std::endl;  streamlog_out(DEBUG3) << "i " ;
    for ( int iii=0 ; iii < 3 ; iii++ ) { streamlog_out(DEBUG3) << DirectionError[iii] << " " ;}  ;streamlog_out(DEBUG3) << std::endl;
 
    FloatVec shape = cluster_->getShape() ; 
	  shape.resize(shape_keys.size());
	  shape[npoints_index]=1.0*n; 
	  shape[sum_wgt_index]=sum_wgt; 
	  shape[sum_wgt2_index]=sum_wgtsqr; 
    shape[sum_wgt4_index]=sum_wgt4; 

    cluster_->setShape( shape ) ;
    cluster_->setPosition(Position);
    cluster_->setPositionError(PositionError);
    cluster_->setITheta(theta);
    cluster_->setIPhi(phi);
    cluster_->setDirectionError(&DirectionError[0]);
 

 
  }


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