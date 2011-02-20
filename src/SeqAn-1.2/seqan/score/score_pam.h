 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: score_pam.h 4615 2009-07-25 07:09:21Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SCORE_PAM_H
#define SEQAN_HEADER_SCORE_PAM_H


//#include <seqan/common.h>
//#include <seqan/score_pam_matrix.h>



//using namespace std;
namespace SEQAN_NAMESPACE_MAIN
{
    
/**
.Spec.Pam:
..cat:Scoring
..summary:Pam scoring matrices.
..general:Class.Score
..signature:Score<TValue, Pam<TSequenceValue, TSpec> >
..param.TValue:Type of the score values.
...default:$int$
..param.TSequenceValue:Type of alphabet underlying the matrix.
...default:$AminoAcid$
..param.TSpec:Origin of starting data underlying the pam matrix computation. The starting data are the mutation probability matrix and the frequencies of occurence of each identifier.
...default:$Pam_Data_Dayhoff_MDM78$
...value:$Pam_Data_Dayhoff_MDM78$, $Pam_Data_Jones_PRI29$, $Pam_Data_Jones_PET91_SWISS15$, $Pam_Data_Jones_PET91_SWISS22$, $Pam_Data_Jones_All_Membrane$, $Pam_Data_Jones_Single_Membrane$, $Pam_Data_Jones_Multi_Membrane$
..remarks.text:This class computes pam matrices for arbitrary alphabets and starting values. The computation precedure includes extrapolation of the MDM (Mutation Probability Matrix) to the desired evolutionary distance in units of 1 PAM, 
computation of the odds scores and the construction of the symmetric and scaled logarithms of odds matrix.	The odds scores thereby represent the chance of a relationship between the observed amino acids at each position versus the chance of coincidental pairing.
If integer scores are produced, each value is rounded separately.
.
.Function.buildPam.param.score.type:Spec.Pam
.Function.score.param.score.type:Spec.Pam
.Function.scoreGapExtend.param.score.type:Spec.Pam
.Function.scoreGapOpen.param.score.type:Spec.Pam
.Function.getScale.param.score.type:Spec.Pam
.Function.getDist.param.score.type:Spec.Pam
.Function.getEntropy.param.score.type:Spec.Pam
.DISABLEDFunction.write.param.data.type:Spec.Pam
.Internal.getPamProperties.param._score.type:Class.Score
.Internal.showPamMatrix.param._score.type:Class.Score
.Internal._setDist.param._score.type:Class.Score
.Internal._setScale.param._score.type:Class.Score
.Internal._getDataPam.param._score.type:Class.Score
.Internal._setEntropy.param._score.type:Class.Score
.Internal._extendAlphabetPam.param._score.type:Class.Score
.Internal._extrapolatePam.param._score.type:Class.Score
..remarks.text:The default invokation creates integer pam matrices for a 250 PAM evolutionary distance from Jones data (1991).
Scoring matrix are stored within the member array "data_pam". Evolutionary distance and entropy are seperatley accessible.
For the usabilty in alignemt alogrithms gap penaltys are stored and externally accessible as well.
..remarks.text: Provided starting data, including frequencies of occurence of each symbol and mutation probabilities between each pair of symbols:
..remarks.text:- Pam_Data_Dayhoff_MDM78: Generated from ATLAS OF PROTEIN SEQUENCE AND STRUCTURE (Margaret Dayhoff et al., 1978)
..remarks.text:- Pam_Data_Jones_PRI29: Generated from PIR Release 29.0 (Jones D.T. et al., 1992)
..remarks.text:- Pam_Data_Jones_PET91_SWISS15: Generated from SWISSPROT Release 15.0 (Jones D.T. et al., 1992)
..remarks.text:- Pam_Data_Jones_PET91_SWISS22: Generated from SWISSPROT Release 22.0 (Jones D.T. et al., 1992)
..remarks.text:- Pam_Data_Jones_All_Membrane: Generated from SWISSPROT Release 23., restricted to transmembranal segments (Jones D.T. et al., 1992)
..remarks.text:- Pam_Data_Jones_Single_Membrane: Generated from SWISSPROT Release 23., restricted to single-spanning transmembranal segments (Jones D.T. et al., 1992)
..remarks.text:- Pam_Data_Jones_Single_Membrane: Generated from SWISSPROT Release 23., restricted to multi-spanning transmembranal segments (Jones D.T. et al., 1992)
*/


template <typename TValue, typename TSequenceValue, typename TSource> 
 class 	Score<TValue, Pam<TSequenceValue, TSource> >		
// TScoureValue: type of source data matrices (e.g. double), TSource: origin of source data (e.g. Dayhoff78)

 {										

public:
	int dist;					// desired PAM distance
	double scaling_factor;		// scale of logarithm of odds matrix
	double entropy;				// information content of final PAM matrix
	double minLength;			// entropy- dependent minimal sequence length for significant scores
	enum{dim = ValueSize<AminoAcid>::VALUE}; // size of alphabet
	TValue data_gap_extend;		// gap opening penalty
	TValue data_gap_open;		// gap extension penalty
	TValue data_pam[dim][dim];	// final PAM matrix


public:

	//////////////////////////////////////////////////////////////////////////////
	// Initialization: either default values or user defined by distance, gap penaltys or scaling factor
	//////////////////////////////////////////////////////////////////////////////

	Score(int _distance = 250, TValue _gap_extend = -1 ):
		data_gap_extend(_gap_extend),
		data_gap_open(_gap_extend)
	{
		buildPam(* this, _distance);
	}
	Score(int _distance, TValue _gap_extend, TValue _gap_open ):
		data_gap_extend(_gap_extend),
		data_gap_open(_gap_open)
	{
		buildPam(* this, _distance);
	}

	/**	
.Memfunc.Score:
..class:Class.Score
..summary:Constructor
..signature:Score(const & other)
..signature:Score(distance, gap [, gap_open])
..param.distance: desired evolutionary distance (int).
...default:250
..param.gap: desired penalty for gap (extension).
...default:-1
..param.gap_open: desired penalty for gap opening.
...default:$gap$
..param.other:Score to be copied
...text:If both gap and gap_open are specified, the total score of a length $n$ gap is $gap_open + (n-1)*gap$.
...note:Usually $gap$, $gap_extend$, and $gap_open$ are negative values.
*/	
		
	Score(Score const & other):
		dist(other.dist),
		scaling_factor(other.scaling_factor),
		entropy (other.entropy),
		data_gap_extend(other.data_gap_extend),
		data_gap_open(other.data_gap_open)
		//data_pam(_getDataPam(* other))
	{
		const TValue * adr_data_pam = _getDataPam(other); 
		//arrayCopyForward(& data_pam[0][0], & data_pam[0][0] + dim * dim, adr_data_pam);
		arrayCopyForward(adr_data_pam, adr_data_pam + dim * dim, & data_pam[0][0]);
		// to do: data_pam kopieren!
	}

		~Score()
	{
	}



/*
	friend inline void
		getPamProperties(Score & _score){
			 int dist = getDist(_score);
			 double scale = getScale(_score);
			 double entropy = getEntropy(_score);
			 cout << "Distance: " << dist << "\n" << "Scaling Factor: " << scale << "\n" << "Entropy: " <<  entropy << "\n\n";	 

		}

	friend inline void
		getPamProperties(Score const & _score){
			 const int dist = getDist(_score);
			 const double scale = getScale(_score);
			 const double entropy = getEntropy(_score);
			 cout << "Distance: " << dist << "\n" << "Scaling Factor: " << scale << "\n" << "Entropy: " <<  entropy << "\n\n";	 

		}

*/
/*DISABLED
.Function.getPamProperties:
..class:Class.Score
..summary:Prints out the properties of the scoring matrix currently contained in "_score". These properties are distance, scaling factor and entropy.
..signature:getPamProperties(& _score)
..signature:getPamProperties(const & _score)
..param._score:Score class instance containing scoring matrix as a member.
*/ 
/*
	friend inline void
		 showPamMatrix(Score & _score){
			 TValue * adr = _getDataPam(_score);
			 for (int i=0;i<dim;i++){
				 for (int j=0; j<dim; j++){
					 cout << *(adr+i*dim+j) << ", ";
				 }
				 cout << "\n";
			 }
		}

 	 friend inline void
		 showPamMatrix(Score const & _score){
			 TValue * adr = _getDataPam(_score);
			 for (int i=0;i<dim;i++){
				 for (int j=0; j<dim; j++){
					 cout << *(adr+i*dim+j) << ", ";
				 }
				 cout << "\n";
			 }
		}
*/

/*DISABLED, use write(.) instead
.Function.showPamMatrix:
..class:Class.Score
..summary:Prints out PAM matrix on display.
..signature:showPamMatrix(& _score)
..signature:showPamMatrix(const & _score)
..param._score:Score class instance containing scoring matrix as a member.
*/
		 
//____________________________________________________________________________



};



	 /**
.Memfunc.~Score:
..class:Class.Score
..summary:Destructor
*/
template <typename TValue, typename TSequenceValue, typename TSource> 
	inline TValue *
		_getDataPam(Score<TValue, Pam<TSequenceValue, TSource> > & _score){
        return & _score.data_pam[0][0];
    }

	
template <typename TValue, typename TSequenceValue, typename TSource> 
	inline const TValue *
		_getDataPam(const Score<TValue, Pam<TSequenceValue, TSource> > & _score){
        return & _score.data_pam[0][0];
    }
	
	/**
.Internal._getDataPam:
..class:Class.Score
..cat:Functions
..summary:Access function returning address of the array containing the scoring matrix.
..signature:_getDataPam(Score & _score)
..param._score:Score class instance containing scoring matrix as a member.
..returns:Address of array containing scoring matrix.
...type:TValue *

*/

/*
	friend inline TValue &
		scoreGapExtend(Score & _score)
	{
		return _score.data_gap_extend;
	}
	friend inline TValue const &
		scoreGapExtend(Score const & _score)
	{
		return _score.data_gap_extend;
	}

	
	friend inline TValue
		scoreGapOpen(Score & _score)
	{
		return _score.data_gap_open;
	}
	friend inline TValue const &
		scoreGapOpen(Score const & _score)
	{
		return _score.data_gap_open;
	}

*/
template <typename TValue, typename TSequenceValue, typename TSource> 
	inline double
		getScale(Score<TValue, Pam<TSequenceValue, TSource> > & _score) {
			return _score.scaling_factor;
		}

template <typename TValue, typename TSequenceValue, typename TSource> 
 	inline double
		getScale(Score<TValue, Pam<TSequenceValue, TSource> > const & _score) {
			return _score.scaling_factor;
		}


/**
.Function.getScale:
..class:Class.Score
..summary:Access function returning the distance dependent $scaling factor$ used for PAM matrix computation
..signature:getScale(score)
..param.score:Score class instance containing scoring matrix as a member.
..returns:Scaling factor used for current PAM matrix computation (double).
*/	
		
template <typename TValue, typename TSequenceValue, typename TSource> 
	inline void
		_setScale(Score<TValue, Pam<TSequenceValue, TSource> > & _score, double _scale) {
			_score.scaling_factor = _scale;
		}


/**
.Internal._setScale:
..class:Class.Score
..cat:Functions
..summary:assigns given value to $scaling_factor$
..signature:_setScale(& _score, _scale)
..param._score:Score class instance in which Pam computation is to be performed.
..param._scale:Is assigned to the member attribute $scaling_factor$ (double).
*/	
		
template <typename TValue, typename TSequenceValue, typename TSource> 
	inline int
		getDist(Score<TValue, Pam<TSequenceValue, TSource> > & _score) {
			return _score.dist;
		}

template <typename TValue, typename TSequenceValue, typename TSource> 
	inline int
		getDist(Score<TValue, Pam<TSequenceValue, TSource> > const & _score) {
			return _score.dist;
		}


/**
.Function.getDist:
..cat:Scoring
..class:Class.Score
..summary:returns the distance on which PAM matrix computation is based.
..signature:getDist(& _score)
..signature:getDist(const & _score)
..param._score:Score class instance containing scoring matrix as a member.
..returns:Distance used for current PAM matrix computation.
*/	
		
template <typename TValue, typename TSequenceValue, typename TSource> 
	inline void
		_setDist(Score<TValue, Pam<TSequenceValue, TSource> > & _score, int _givenDist) {
			_score.dist = _givenDist;
		}

/**
.Internal._setDist:
..class:Class.Score
..cat:Functions
..summary:assigns given value to $dist$
..signature:_setDist(& _score, _givenDist)
..param._score:Score class instance in which Pam computation is to be performed.
..param._givenDist:Is assigned to the member attribute $dist$.
*/
template <typename TValue, typename TSequenceValue, typename TSource> 
	inline double
		getEntropy(Score<TValue, Pam<TSequenceValue, TSource> > & _score) {
			return _score.entropy;
		}


template <typename TValue, typename TSequenceValue, typename TSource> 
	inline double
		getEntropy(Score<TValue, Pam<TSequenceValue, TSource> > const & _score) {
			return _score.entropy;
		}

/**
.Function.getEntropy:
..cat:Scoring
..class:Class.Score
..summary:returns entropy of the PAM matrix.
..signature:getEntropy(& _score)
..signature:getEntropy(const & _score)
..param._score:Score class instance containing scoring matrix as a member.
..returns:Entropy of the matrix currently hold in Score $_score$ (double).
*/
template <typename TValue, typename TSequenceValue, typename TSource> 
	inline void
		_setEntropy(Score<TValue, Pam<TSequenceValue, TSource> > & _score, double H) {
			_score.entropy = H;
		}
/**
.Internal._setEntropy:
..class:Class.Score
..cat:Functions
..summary:assigns given value to $entropy$.
..signature:_setEntropy(Score & _score, double H)
..param._score:Score class instance containing scoring matrix as a member.
..param.H:Is assigned to the member attribute $entropy$ (double).
*/


//////////////////////////////////////////////////////////////////////////////
// Define some types for easier instantiation of PAM matrix
//////////////////////////////////////////////////////////////////////////////

typedef Score <int, Pam<AminoAcid, Pam_Data_Dayhoff_MDM78> >	PamDayhoff;
typedef Score <int, Pam<AminoAcid, Pam_Data_Jones_PRI29> >	PamJones;
typedef Score <int, Pam<AminoAcid, Pam_Data_Jones_PET91_SWISS15> > PamJonesSWISS15;
typedef Score <int, Pam<AminoAcid, Pam_Data_Jones_PET91_SWISS22> > PamJonesSWISS22;
typedef Score <int, Pam<AminoAcid, Pam_Data_Jones_All_Membrane> > PamMembrane;
typedef Score <int, Pam<AminoAcid, Pam_Data_Jones_Single_Membrane> > PamSingleMembrane;
typedef Score <int, Pam<AminoAcid, Pam_Data_Jones_Multi_Membrane> > PamMultiMembrane;





//////////////////////////////////////////////////////////////////////////////
//					 Functions for PAM matrix computation								 
//////////////////////////////////////////////////////////////////////////////
	template <typename TValue, typename TSequenceValue, typename TSource>
void buildPam(Score<TValue, Pam<TSequenceValue, TSource> > & _score, int _givenDist)
{
/**
.Function.buildPam:
..cat:Scoring
..class:Class.Score
..summary:Coordinates computation by invokating extrapolation, alphabet extension, rounding, scaling, type conversion and symmetrization.
..signature:buildPam(_score, _givenDist)
..param._score:Score class instance in which PAM computation is to be performed. 
..param._givenDist:Evolutionary distance underlying PAM computation
*/
	SEQAN_CHECKPOINT
	_TempMembersPam<TValue, TSequenceValue> my_members;
	my_members.min_score = TValue(10e5);

	_setDist(_score, _givenDist);
	_computePamScale(_score);
	_starting_data_pam(my_members, _score);
	
	_extrapolatePam(my_members, _score);
	_extendAlphabetPam(_score, my_members);
	_computeLogOddsPam(my_members);
	_finishPam(_score, my_members);
	_computeEntropyPam(_score, my_members);
//____________________________________________________________________________

}



template <typename TValue, typename TSource, typename TSequenceValue>
void _extrapolatePam (_TempMembersPam<TValue, TSequenceValue> & _member, Score<TValue, Pam<TSequenceValue, TSource> > _score){

/**
.Internal._extrapolatePam:
..class:Class.Score
..cat:Functions
..summary:Extrapolates the mutation probability matrix to the desired evolutionary distance. The extrapolation is performed parsimoniously by determining the least number of squares necessary to reach the desired distance.
..signature:(_member, _score)
..param._member:Temporary struct containing starting data and intermediate matrix products as members. 
...type:_TempMembersPam
..param._score:Score class instance in which PAM computation is to be performed. 
*/	
	
	SEQAN_CHECKPOINT
	// extrapolate matrix to distance between 2 and 512
	enum 
	{
		dim = ValueSize<TSequenceValue>::VALUE
	};
	int _dist = getDist(_score);
	 SEQAN_TASSERT2( (_dist>2) && (_dist < 512), "Distance must be between 2 and 512 for protein comparisons!")
	
	
	double temp_prob[dim][dim];		// temporary target for intermediate results
	double temp_product[dim][dim];	// help array for procedure of matrix multiplication
	int bit_count = (sizeof(int)* 8);	// Anzahl der Bitvergleiche

	for (int i=0; i<dim; i++)	// temporary target is initially an identity-matrix
	{		
		for (int j=0; j<dim; j++) 
		{
			if (i==j) {	temp_prob[i][j] = 1.;}
			else {temp_prob[i][j] = 0.;}
		}
	}	
	for (int i=1; i<=bit_count; i++) {
		if (_dist != 0) // compare values and stop when distance is shifted to zero
		{	
			if (_dist & 1) // compare bitwise
			{		
				//_multiply_self(temp_prob,_extrapolateMe);
				for (int i=0; i<dim; i++){
					for (int j=0; j<dim; j++){
						temp_product[i][j] = 0.;
						for (int k=0; k<dim; k++) {
							temp_product[i][j] += temp_prob[i][k] * _member.trans_prob[k][j];
						}
					}
				}
				arrayCopyForward(* temp_product,* temp_product+dim*dim,* temp_prob);

			}
			for (int i=0; i<dim; i++){
					for (int j=0; j<dim; j++){
						temp_product[i][j] = 0.;
						for (int k=0; k<dim; k++) {
							temp_product[i][j] += _member.trans_prob[i][k] * _member.trans_prob[k][j];
						}
					}
				}
			arrayCopyForward(* temp_product,* temp_product+dim*dim,* _member.trans_prob);

			_dist = _dist >> 1;
		}	
	}

	for (int i=0; i<dim-4; i++){
		for (int j=0; j<dim-4; j++) {
			temp_product[i][j] = temp_product[j][i]	= ((temp_product[i][j] + temp_product[j][i])/2);

		}
	}
	arrayCopyForward(* temp_prob,* temp_prob+dim*dim,* _member.trans_prob);
//____________________________________________________________________________
}









template <typename TValue, typename TSource>
void _extendAlphabetPam(Score<TValue, Pam<AminoAcid, TSource> > &, _TempMembersPam<TValue,AminoAcid> & _member){
// extension of amino acid alphabet: known location for B and Z 
	
/**
.Internal._extendAlphabetPam:
..class:Class.Score
..cat:Functions
..summary:If the AminoAcid alphabet is used, this function extends the matrix by the unkown identifiers 'B', 'Z', '*', 'X'.
..signature:_extendAlphabetPam(& _score, _member)
..param._score:Score class instance in which PAM computation is to be performed. 
..param._member:Temporary struct containing starting data and intermediate matrix products as members. 
...type:_TempMembersPam
*/	
	
	const int dim = _member.dim;

	_member.single_prob[dim-4] = _member.single_prob[2] + _member.single_prob[3];
	_member.single_prob[dim-3] = _member.single_prob[5] + _member.single_prob[6];

	for (int i=0; i<dim-2; i++) {
		_member.trans_prob[i][dim-4] = ((_member.single_prob[2]*_member.trans_prob[i][2])+(_member.single_prob[3]*_member.trans_prob[i][3]))/(_member.single_prob[dim-4]) ;
		_member.trans_prob[i][dim-3] = ((_member.single_prob[5]*_member.trans_prob[i][5])+(_member.single_prob[6]*_member.trans_prob[i][6]))/(_member.single_prob[dim-3]) ;
		_member.trans_prob[dim-4][i] = _member.trans_prob[2][i] + _member.trans_prob[3][i];
		_member.trans_prob[dim-3][i] = _member.trans_prob[5][i] + _member.trans_prob[6][i];

	}
//____________________________________________________________________________

}

template <typename TValue, typename TSequenceValue, typename TSource>
void _extendAlphabetPam(Score<TValue, Pam<TSequenceValue, TSource> > & _score, _TempMembersPam<TValue, TSequenceValue> & _member){
// Dummy for arbitrary alphabets --> no extension
//____________________________________________________________________________

}











template <typename TValue, typename TSequenceValue>
void _computeLogOddsPam(_TempMembersPam<TValue, TSequenceValue> & _member){
/**
.Internal._computeLogOddsPam:
..class:Class.Score
..cat:Functions
..summary:Computes logarithms of odds score of the extrapolated mutation probability matrix. Each odds score represents the chance of evolutionary relationship at the given distance versus the chance of no such relationship.
..signature:_computeLogOddsPam(_member)
..param._member:Temporary struct containing starting data and additional information as members. 
...type:_TempMembersPam
*/
	SEQAN_CHECKPOINT

	// compute natural (base e) logarithm of odds 
	const int dim = _member.dim;
	double log_numerator = 1.; double log_denominator = 1.;
	for (int i=0; i<dim-2 ; i++){
			for (int j=0; j<dim-2 ; j++){
				log_numerator = std::log(_member.trans_prob[i][j]); 
				log_denominator = std::log(_member.single_prob[i]);
				_member.log_prob[i][j] = (log_numerator - log_denominator);
			}
		}
//____________________________________________________________________________
}


template <typename TValue, typename TSequenceValue, typename TSource>
void _finishPam(Score<TValue, Pam<TSequenceValue, TSource> > & _score, _TempMembersPam<TValue, TSequenceValue> & _member){
// Derive Pam scores for arbitrary alphabets. Only the unextended alphabet is used.

	const int dim = _member.dim;
	int terminator = dim-1, xIndex = dim-2;
	double whole_sym_entry = 0.;
	int int_entry = 0;
	double _scale = getScale(_score);
	TValue * adr_data_pam = _getDataPam(_score); 

	//////////////////////////////////////////////////////////////////////////////
	// fill matrix of all identifiers: rounding, scaling and type conversion
	//////////////////////////////////////////////////////////////////////////////
	for (int i=0; i<dim-2; i++) {
			for (int j=0; j<=i; j++) {
				if (i != j) { // nondiagonal elements: symmetry adjustment and scaling
					whole_sym_entry = ((_member.log_prob[i][j]+_member.log_prob[j][i])*_scale)/2.;
				}
				else { // diagonal elements: just scaling
					whole_sym_entry = _member.log_prob[i][j]*_scale;	
				}
				int_entry = roundConvert(whole_sym_entry, _score);	// round correctly
				* (adr_data_pam + i*dim + j) = 	*(adr_data_pam + j*dim + i) = int_entry;
				if (int_entry < _member.min_score) {_member.min_score = int_entry;	}
			}
		}
	//____________________________________________________________________________


	//////////////////////////////////////////////////////////////////////////////
	// Compute score for 'unkown' character 'X'
	//////////////////////////////////////////////////////////////////////////////

	// Score for X versus standard amino acids
	double xScoreSelf = 0;		// Score for X versus X
	double xScoreChange = 0;	// Score for X versus any standard amino acid

	for (int i=0; i<dim-2; i++){
		xScoreChange = 0;
		for (int j=0; j<dim-4; j++){
			// Compute Score for X versus standard amino acid
			xScoreChange += _member.single_prob[j] * (*(adr_data_pam + i*dim + j));
			// Compute Score for X versus X
			xScoreSelf += _member.single_prob[j] * _member.single_prob[i] * (*(adr_data_pam + i*dim + j));
		}

		// Assign score for X versus standard amino acid to matrix
		*(adr_data_pam + xIndex*dim + i)   =   *(adr_data_pam + i*dim + xIndex) = roundConvert(xScoreChange, _score);
	}

	// Assign score for X versus X to matrix
	*(adr_data_pam + xIndex*dim + xIndex) = roundConvert(xScoreSelf, _score);


	//____________________________________________________________________________




	//////////////////////////////////////////////////////////////////////////////
	// Compute score for 'terminator' character '*'
	//////////////////////////////////////////////////////////////////////////////

	for (int i=0; i<dim; i++) {
  			*(adr_data_pam + terminator*dim + i )  = *(adr_data_pam + i*dim + terminator )	=   _member.min_score;
		}
	*(adr_data_pam + terminator*dim + terminator ) = 1;
	//____________________________________________________________________________

}



template <typename TValue, typename TSource>
void _finishPam(Score<TValue, Pam<AminoAcid, TSource> > & _score, _TempMembersPam<TValue, AminoAcid> & _member){
// Derive Pam scores for the AminoAcid alphabet. Computation of the 'unkown' character 'X' is extended to the additional identifiers 'B' and 'Z'.

/**
.Internal._finishPam:
..class:Class.Score
..cat:Functions
..summary:Creates final Pam matrix from log odds matrix. Each value is scaled by multiplication with $scaling_factor$. If necessary, the score values are rounded before inserted into the final score matrix. In case of the AminoAcid alphabet, values for the IUPAC identifiers 'B', 'Z', 'X', '*' are derived.
..signature:_finishPam(& _score, & _member)
..param._score:Score class instance in which PAM computation is to be performed. 
..param._member:Temporary struct containing starting data and additional information as members. 
...type:_TempMembersPam
*/
	const int dim = _member.dim;
	int terminator = dim-1, xIndex = dim-2;
	double whole_sym_entry = 0.;
	TValue int_entry = 0;
	double _scale = getScale(_score);
	TValue * adr_data_pam = _getDataPam(_score); 

	//////////////////////////////////////////////////////////////////////////////
	// fill matrix for all 22 IUPAC amino acids: rounding, scaling and type conversion
	//////////////////////////////////////////////////////////////////////////////
	for (int i=0; i<xIndex; i++) {
			for (int j=0; j<=i; j++) {
				if (i != j) { // nondiagonal elements: symmetry adjustment and scaling
					whole_sym_entry = ((_member.log_prob[i][j]+_member.log_prob[j][i])*_scale)/2.;
				}
				else { // diagonal elements: just scaling
					whole_sym_entry = _member.log_prob[i][j]*_scale;	
				}
				int_entry = roundConvert(whole_sym_entry, _score);	// round correctly
				* (adr_data_pam + i*dim + j) = 	*(adr_data_pam + j*dim + i) = int_entry;
				if (int_entry < _member.min_score) {_member.min_score = int_entry;	}
			}
		}
	//____________________________________________________________________________


	//////////////////////////////////////////////////////////////////////////////
	// Compute score for 'unkown' character 'X'
	//////////////////////////////////////////////////////////////////////////////

	// Score for X versus standard amino acids
	double xScoreSelf = 0;		// Score for X versus X
	double xScoreChange = 0;	// Score for X versus any standard amino acid
	double xbScore = 0.;		// Score for X versus B 
	double xzScore = 0.;		// Score for X versus Z

	for (int i=0; i<dim-4; i++){
		xScoreChange = 0;
		for (int j=0; j<dim-4; j++){
			// Compute Score for X versus standard amino acid
			xScoreChange += _member.single_prob[j] * (*(adr_data_pam + i*dim + j));
			// Compute Score for X versus X
			xScoreSelf += _member.single_prob[j] * _member.single_prob[i] * (*(adr_data_pam + i*dim + j));
		}
		// Compute Score for X versus B
		xbScore += (_member.single_prob[2]*_member.single_prob[i]*(*(adr_data_pam + 2*dim + i)) +  _member.single_prob[3]*_member.single_prob[i]*(*(adr_data_pam + 3*dim + i)))/_member.single_prob[dim-4];
	   	// Compute Score for X versus Z
		xzScore += (_member.single_prob[5]*_member.single_prob[i]*(*(adr_data_pam + 5*dim + i)) +  _member.single_prob[6]*_member.single_prob[i]*(*(adr_data_pam + 6*dim + i)))/_member.single_prob[dim-3];

		// Assign score for X versus standard amino acid to matrix
		*(adr_data_pam + xIndex*dim + i)   =   *(adr_data_pam + i*dim + xIndex) = roundConvert(xScoreChange, _score);
	}

	// Assign score for X versus X to matrix
	*(adr_data_pam + xIndex*dim + xIndex) = roundConvert(xScoreSelf, _score);
	// Assign score for X versus B to matrix
	*(adr_data_pam + (dim-4)*dim + xIndex) = *(adr_data_pam + xIndex*dim + (dim-4)) = roundConvert(xbScore, _score);
	// Assign score for X versus Z to matrix
	*(adr_data_pam + (dim-3)*dim + xIndex) = *(adr_data_pam + xIndex*dim + (dim-3)) = roundConvert(xzScore, _score);


	//____________________________________________________________________________




	//////////////////////////////////////////////////////////////////////////////
	// Compute score for 'terminator' character '*'
	//////////////////////////////////////////////////////////////////////////////

	for (int i=0; i<dim; i++) {
  			*(adr_data_pam + terminator*dim + i )  = *(adr_data_pam + i*dim + terminator )	=   _member.min_score;
		}
	*(adr_data_pam + terminator*dim + terminator ) = 1;
	//____________________________________________________________________________

}

template <typename TValue, typename TSequenceValue, typename TSource>
void _computeEntropyPam(Score<TValue, Pam<TSequenceValue, TSource> > & _score , _TempMembersPam<TValue, TSequenceValue> & _member ){
	//////////////////////////////////////////////////////////////////////////////
	// General procedure for arbitrary alphabets: use whole range of matrix values for entropy computation
	//////////////////////////////////////////////////////////////////////////////
	const int dim = _member.dim;
	double H = 0.;
	for (int i=0;i<dim;i++) {
		for (int j=0;j<dim;j++) {
			H += _member.single_prob[i]* _member.trans_prob[i][j] * _member.log_prob[i][j];		}
	}
	// Divide entropy by log(2) to give information in units of bits:
	_setEntropy(_score, H/std::log(2.));
	//____________________________________________________________________________

}





template <typename TValue, typename TSource>
void _computeEntropyPam(Score<TValue, Pam<AminoAcid, TSource> > & _score , _TempMembersPam<TValue, AminoAcid> & _member ){
	//////////////////////////////////////////////////////////////////////////////	
	// Specialzation for IUPAC matrices: only 20 single (uncombined) amino acid scores usable for entropy computation
	//////////////////////////////////////////////////////////////////////////////

/**
.Internal._computeEntropyPam:
..class:Class.Score
..cat:Functions
..summary:Computes entropy of current Pam matrix.
..signature:(& _score , & _member ){
..param._score:Score class instance containing scoring matrix as a member.
..param._member:Temporary struct containing starting data and additional information as members. 
...type:_TempMembersPam
*/


	const int dim = _member.dim;
	double H = 0.;
	for (int i=0;i<dim-4;i++) {	// only iterate through 20x20 matrix of single amino acid scores
		for (int j=0;j<dim-4;j++) {
			H += _member.single_prob[i]* _member.trans_prob[i][j] * _member.log_prob[i][j];
		}
	}
	// Divide entropy by log(2) to give information in units of bits:
	_setEntropy(_score, H/std::log(2.));	

	//____________________________________________________________________________

}




template <typename TValue, typename TSequenceValue, typename TSource>
void _computePamScale(Score<TValue, Pam<TSequenceValue, TSource> > & _score) {
	//////////////////////////////////////////////////////////////////////////////		
	// Scaling factor for creation of final score. 
	// - Numerator: enlarge weight on each score with increasing distance to prevent information loss
	// - Denominator: divide natural (base e) log odds score by natural logarithm of 2 
	//   to create log odds scores to base 2	(interpretable as bit units)
	//////////////////////////////////////////////////////////////////////////////		
	SEQAN_CHECKPOINT

	double 	numerator = 1.;
	int dist = getDist(_score);
	if (dist >= 470) {numerator = 7.;}
	else if (dist >= 409) {numerator = 6.;}
	else if (dist >= 341) {numerator = 5.;}
	else if (dist >= 262) {numerator = 4.;}
	else if (dist >= 170) {numerator = 3.;}
	else if (dist >= 1) {numerator = 2.;}
	_setScale(_score, numerator/std::log(2.));

	//____________________________________________________________________________

}




//////////////////////////////////////////////////////////////////////////////
// Help function for type conversion: only round when integer matrix desired
//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TValue, typename TSequenceValue >
TValue roundConvert(TSource to_round, Score<TValue, TSequenceValue> _roundType){
	return TValue(to_round);

}

template <typename TSource, typename TSequenceValue>
int roundConvert(TSource to_round, Score<int, TSequenceValue>){
	return int(floor(to_round + 0.5));
}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



template <typename TFile, typename TValue, typename TSequenceValue, typename TSpec, typename TMeta>
void
write(TFile & fl,
	  Score<TValue, Pam<TSequenceValue, TSpec> > & sc,
	  TMeta & meta)
{
	_writeScoringMatrix<TSequenceValue>(fl, _getDataPam(sc), meta);
}

template <typename TFile, typename TValue, typename TSequenceValue, typename TSpec>
void
write(TFile & fl,
	  Score<TValue, Pam<TSequenceValue, TSpec> > & sc)
{
	char meta[100];
	sprintf(meta, " PAM %i substitution matrix\n\n Scaling: %f\n Entropy: %f\n\n", getDist(sc), getScale(sc), getEntropy(sc));
	write(fl, sc, meta);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSequenceValue, typename TSource, typename TVal1, typename TVal2> 
inline TValue &
score(Score<TValue, Pam<TSequenceValue, TSource> > & _score, 
	  TVal1 as_1, 
	  TVal2 as_2)
{
	TValue * data_pam = _getDataPam(_score);
	TValue & pairscore = *(data_pam+(ordValue(as_1)*_score.dim)+ordValue(as_2));
	return pairscore;

}
template <typename TValue, typename TSequenceValue, typename TSource, typename TVal1, typename TVal2> 
inline TValue const &
score(Score<TValue, Pam<TSequenceValue, TSource> > const & _score, 
	  TVal1 as_1, 
	  TVal2 as_2)
{
	const TValue * data_pam = _getDataPam(_score);
	//return _score.data_pam[as_1][as_2];
	const TValue & pairscore = *(data_pam+(ordValue(as_1)*_score.dim)+ordValue(as_2));
	return pairscore;
}


//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
