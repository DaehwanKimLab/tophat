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
  $Id: markov_model.h 2009-03-10 utro@math.unipa.it srombo@deis.unical.it $
 ==========================================================================*/

//////////////////////////////////////////////////////////////////////////////

#ifndef SEQAN_HEADER_STATISTICAL_INDEX_MARKOV_MODEL_H
#define SEQAN_HEADER_STATISTICAL_INDEX_MARKOV_MODEL_H

namespace seqan
{

unsigned int _getStringIndexFromMatrix(unsigned int ri, unsigned int ci, unsigned int ncol);

/**
.Class.MarkovModel:
..summary:Gives a suitable representation of a Marcov Chain.
..cat:MarkovModel
..signature:MarkovModel<TAlphabet,[TFloat, TSpec]>
..param.TAlphabet:The alphabet type
..param.TFloat:The type of the exploited arrays
..param.TSpec:The MarkovModel type
.Memvar.MarkovModel#order:
..class:Class.MarkovModel
..summary:The MarkovModel order
...type:nolink:int
.Memvar.MarkovModel#transition:
..class:Class.MarkovModel
..summary:The transition matrix
...type:Class.String
.Memvar.MarkovModel#stationaryDistribution:
..class:Class.MarkovModel
..summary:The vector of character distribution
...type:Class.String
.Memfunc.MarkovModel#MarkovModel:
..class:Class.MarkovModel
..summary:Constructor
..signature:MarkovModel(order_)
..param.order_:The order of the MarkovModel.
.Memfunc.MarkovModel#build:
..class:Class.MarkovModel
..summary:Given a training set, computes the transition matrix, the character stationary distributions and the auxiliary information that give raise to an instance of MarkovModel
..signature:build(strings)
..param.strings:The training set.
...type:Class.StringSet

.Memfunc.MarkovModel#set:
..class:Class.MarkovModel
..summary: Given e transition matrix, sets it as transition matrix of the MarkovModel and computes (if it is not available) the  vector of character distributions and the auxiliary information
..signature:set(iTransition)
..set(iTransition, iStationaryDistribution)
..param.iTransition:The transition matrix.
...type:Class.String
..param.iStationaryDistribution:The vector of character distributions.
...type:Class.String

.Memfunc.MarkovModel#emittedProbability:
..class:Class.MarkovModel
..summary:Computes the probability that a string (or a set of strings) is emitted by the MarkovModel.
..signature:emittedProbability(string)
..signature:emittedProbability(stringSet)
..param.string:The string whose emission probability has to be computed.
...type:Class.String
..param.stringSet:The set of strings whose emission probability has to be computed.
...type:Class.StringSet
..returns:A TFloat representing the emission probability.

.Memfunc.MarkovModel#write:
..class:Class.MarkovModel
..summary: Stores an instance of MarkovModel on a file
..signature:write(file)
..param.file:The file on which storing the MarkovModel.

.Memfunc.MarkovModel#read:
..class:Class.MarkovModel
..summary: Loads an instance of MarkovModel from a file
..signature:read(file)
..param.file:The file from which loading the MarkovModel.
*/
template <typename TAlphabet, typename TFloat = double, typename TSpec = Default>
class MarkovModel
{

public:
	typedef String<TFloat> TMatrix;
	typedef String<TFloat> TVector;

	unsigned int order;
	TMatrix transition;
	TVector stationaryDistribution;
	//The following matrices are only for internal use of the class
	TMatrix _q;
	TMatrix _qppp;
	TMatrix _qppqpp;


	MarkovModel(unsigned int order_) :
		order(order_)
	{
		unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
		unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);
		unsigned int const matrix_size = column_size * column_size;

		clear(transition);
		clear(stationaryDistribution);
		fill(transition, matrix_size, 0);
		resize(stationaryDistribution, column_size);
	}


	///////////////////////////////////////////////////////////////
	///// BUILD THE MODEL
	///////////////////////////////////////////////////////////////

	template <typename TStringSet>
	void build(TStringSet & strings)
	{
		typedef Index<TStringSet, Index_QGram<SimpleShape> > TIndex;
		typedef typename Fibre<TIndex, QGram_Dir>::Type TDir;
		typedef typename Iterator<TDir, Standard>::Type TIter;
		unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
		unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);

		TIndex ind(strings);
		resize(indexShape(ind), order + 1);
		indexRequire(ind, QGram_SADir());

		TIter itBegin = begin(indexDir(ind), Standard());
		TIter itEnd = end(indexDir(ind), Standard()) - 1;

		String<TAlphabet> qgram;
		Shape<TAlphabet, SimpleShape> orderShape;
		resize(orderShape, order);
		for(TIter itA = itBegin; itA != itEnd; ++itA) {
			unhash(qgram, itA - itBegin, weight(indexShape(ind)));
			value(transition, hash(orderShape, begin(qgram)) * column_size + hash(orderShape, begin(qgram) + 1)) = *(itA+1) - *itA;
		}

		//normalization
		for(unsigned int row = 0; row < column_size;++row) {
		TFloat sum = 0;
		for(unsigned int col = 0; col < column_size;++col) {
			sum += value(transition, (row * column_size) + col);
		}
		if (sum != 0) {
				for(unsigned int col = 0; col < column_size;++col) {
					value(transition, row * column_size + col) /= sum;
				}
			}
		}

		String<TFloat> temp = transition;
		//initialize a variable t representing a good trashold to extimate the vector
		//after multiplying e times the transition matrix with itself
		unsigned int t=6;
		for (unsigned int i=0; i<t; i++){
			temp=_matricialProduct(temp,temp,column_size);
		}

		for (unsigned int i=0; i<column_size; i++){
			value(stationaryDistribution,i)=value(temp, i);
		}

		_computeAuxiliaryMatrices();

	}


	///////////////////////////////////////////////////////////////
	///// EMITTEDPROBABILITY
	///////////////////////////////////////////////////////////////


	template <typename TStringSet>
	TFloat emittedProbability(TStringSet &string)
	{
		TFloat p = 0;

		for(unsigned int i=0; i<length(string); i++)
		{
			p+= emittedProbability(getValueById(string, i));
		}

		return p;
	}


    TFloat emittedProbability(String<TAlphabet> &string)
	{
		Shape<TAlphabet, SimpleShape> orderShape;
		resize(orderShape, order);
		unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
		unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);

		int row = hash(orderShape,begin(string));
		TFloat p = value(stationaryDistribution,row);

		for(unsigned int i=1; i<length(string)-order+1; i++)
		{
			int column=hash(orderShape,begin(string)+i);
			p*=value(transition,_getStringIndexFromMatrix(row,column,column_size));
			row = column;
		}

		return p;
	}

	///////////////////////////////////////////////////////////////
	///// SET THE MODEL
	///////////////////////////////////////////////////////////////


	void set(String<TFloat> &iTransition)
	{
		unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
		unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);

		transition = iTransition;

		String<TFloat> temp = transition;
		//initialize a variable t representing a good trashold to extimate the vector
		//after multiplying e times the transition matrix with itself
		unsigned int t=6;
		for (unsigned int i=0; i<t; i++){
			temp=_matricialProduct(temp,temp,column_size);
		}

		for (unsigned int i=0; i<column_size; i++){
			value(stationaryDistribution,i)=value(temp, i);
		}

		_computeAuxiliaryMatrices();
	}



	void set(String<TFloat> &iTransition, String<TFloat> &iStationaryDistribution)
	{
		transition = iTransition;

		stationaryDistribution = iStationaryDistribution;

		_computeAuxiliaryMatrices();
	}


	///////////////////////////////////////////////////////////////
	///// WRITE
	///////////////////////////////////////////////////////////////


	void write(FILE *file)
	{
		unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
		unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);

		//write the transition matrix
		for(unsigned int row=0; row<column_size; row++){
			for(unsigned int col=0; col<column_size; col++){
			  fprintf(file,"%f ",value(transition, _getStringIndexFromMatrix(row,col,column_size)));
			}
			fprintf(file,"\n");
		}
		//write the stationary distribution vector
		for(unsigned int row=0; row<column_size; row++){
			  fprintf(file,"%f ",value(stationaryDistribution, row));
			}
		fprintf(file,"\n");

		if(length(_q)){
			//write the auxiliary matrix
			for(unsigned int row=0; row<column_size; row++){
				for(unsigned int col=0; col<column_size; col++){
					fprintf(file,"%f ",value(_q, _getStringIndexFromMatrix(row,col,column_size)));
				}
				fprintf(file,"\n");
			}

			for(unsigned int row=0; row<column_size; row++){
				for(unsigned int col=0; col<column_size; col++){
				  fprintf(file,"%f ",value(_qppp, _getStringIndexFromMatrix(row,col,column_size)));
				}
				fprintf(file,"\n");
			}

			for(unsigned int row=0; row<column_size; row++){
				for(unsigned int col=0; col<column_size; col++){
				fprintf(file,"%f ",value(_qppqpp, _getStringIndexFromMatrix(row,col,column_size)));
				}
				fprintf(file,"\n");
			}
		}
	}


	///////////////////////////////////////////////////////////////
	///// READ
	///////////////////////////////////////////////////////////////

	void read(FILE *file)
	{
		unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
		unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);
		unsigned int const matrix_size = column_size * column_size;

		//read the transition matrix
		for(unsigned int row=0; row<column_size; row++){
			for(unsigned int col=0; col<column_size; col++){
			  fscanf(file,"%lf ", & value(transition, _getStringIndexFromMatrix(row,col,column_size)));
			}
			fscanf(file,"\n");
		}
		//read the stationary distribution vector
		for(unsigned int row=0; row<column_size; row++){
			  fscanf(file,"%lf ",&value(stationaryDistribution, row));
			}
		fscanf(file,"\n");

		if (!feof(file)){

			resize(_q, matrix_size);

			//read the auxiliary matrix
			for(unsigned int row=0; row<column_size; row++){
				for(unsigned int col=0; col<column_size; col++){
					fscanf(file,"%lf ",&value(_q, _getStringIndexFromMatrix(row,col,column_size)));
				}
				fscanf(file,"\n");
			}

			resize(_qppp, matrix_size);
			for(unsigned int row=0; row<column_size; row++){
				for(unsigned int col=0; col<column_size; col++){
				  fscanf(file,"%lf ",&value(_qppp, _getStringIndexFromMatrix(row,col,column_size)));
				}
				fscanf(file,"\n");
			}

			resize(_qppqpp, matrix_size);
			for(unsigned int row=0; row<column_size; row++){
				for(unsigned int col=0; col<column_size; col++){
					fscanf(file,"%lf ",&value(_qppqpp, _getStringIndexFromMatrix(row,col,column_size)));
				}
				fscanf(file,"\n");
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////
	///// COMPUTE THE AUXILIARY MATRICES FOR THE VARIANCE AND Z-SCORE COMPUTATION
	/////////////////////////////////////////////////////////////////////////////

	/*
		.Memfunc.MarkovModel#_computeAuxiliaryMatrices:
		..class:Class.MarkovModel
		..summary:Computes the auxiliary information for statistical indices computation
		..signature:_computeAuxiliaryMatrices()
	*/

	void _computeAuxiliaryMatrices()
	{

		clear(_q);
		clear(_qppp);
		clear(_qppqpp);

		unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
		unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);
		unsigned int const matrix_size = column_size * column_size;
		TMatrix I;
		TMatrix Ip;
		fill(Ip,matrix_size,0);
		fill(I,matrix_size,0);

		for(unsigned int i=0; i<column_size; i++){
			value(I,_getStringIndexFromMatrix(i,i,column_size))=1;
			 for (unsigned int j=0; j<column_size; j++)
			    value(Ip,_getStringIndexFromMatrix(i,j,column_size))=value(stationaryDistribution,j);
		}



		_q=_computeInverseMatrix(_matricialSum(_matricialDifference(transition, I, column_size), Ip, column_size), column_size);

		_qppp=_matricialProduct(_q,transition,column_size);
		_qppqpp = _matricialProduct(_matricialProduct(_qppp,transition, column_size),_q,column_size);
		for(unsigned int i=1; i<order; i++){
			_qppp=_matricialProduct(_qppp,transition,column_size);
			_qppqpp=_matricialProduct(_qppqpp,transition,column_size);
		}

	}

};


//CODE TO MANAGE MATRICES

/*
.Function._getMatrixRowFromString:
..summary:Gets the row index of a matrix stored as a string, given a string position
..signature:_getMatrixRowFromString(stringPosition, ncol)
..param.stringPosition:The string position.
...type:nolink:unsigned int
..param.ncol:The number of columns of the matrix.
...type:nolink:unsigned int
..returns:The matrix row index.
*/

 unsigned int _getMatrixRowFromString(unsigned int stringPosition, unsigned int ncol)
{
	return (int)(stringPosition/ncol);
}


/*
.Function._getMatrixColumnFromString:
..summary:Gets the column index of a matrix stored as a string, given a string position
..signature:_getMatrixRowFromString(stringPosition, ncol)
..param.stringPosition:The string position.
...type:nolink:unsigned int
..param.ncol:The number of columns of the matrix.
...type:nolink:unsigned int
..returns:The matrix row index.
*/

 unsigned int _getMatrixColumnFromString(unsigned int stringPosition, unsigned int ncol)
{
	unsigned int row_index=_getMatrixRowFromString(stringPosition, ncol);
	return (int)(stringPosition-(ncol*row_index));
}


/*
.Function._getStringIndexFromMatrix:
..summary:Gets the string index of a matrix stored as a string, given the row and column indexes
..signature:_getStringIndexFromMatrix(ri, ci, ncol)
..param.ri:The row index.
...type:nolink:unsigned int
..param.ci:The column index.
...type:nolink:unsigned int
..param.ncol:The number of columns of the matrix.
...type:nolink:unsigned int
..returns:The string index of the matrix.
*/

 unsigned int _getStringIndexFromMatrix(unsigned int ri, unsigned int ci, unsigned int ncol)
{
	return ((ri*ncol)+ci);
}


/*
.Function._matricialProduct:
..summary:Computes the matricial product between two matrixes
..signature:_matricialProduct(matrix1,matrix2,n)
..param.matrix1:The first matrix.
...type:String<TAlphabet, TSpec1>&
..param.matrix2:The second matrix.
...type:String<TAlphabet, TSpec1>&
..param.n:The number of rows and columns of the first and of the second matrix, resp..
...type:nolink:unsigned int
..returns:The products of the two matrices (another matrix).
..remarks:The number of rows of matrix1 must be equal to the number of columns of matrix2.
*/

template <typename TAlphabet, typename TSpec1, typename TSpec2>
String<TAlphabet> _matricialProduct(String<TAlphabet, TSpec1>& matrix1,
					  String<TAlphabet, TSpec2>& matrix2, unsigned int n)
{
	//SEQAN_ASSERT(((length(matrix1)%n)==0)&&((length(matrix2)%n)==0));
	unsigned int nrow1=length(matrix1)/n;
	unsigned int ncol2=length(matrix2)/n;
	String<TAlphabet> result;
	unsigned int lenghtnew=nrow1*ncol2;
	fill(result,lenghtnew,0);

	for(unsigned int row = 0; row < nrow1; row++){
			for(unsigned int col = 0; col < ncol2; col++){
				for(unsigned int colRes = 0; colRes < n; colRes++){

					value(result, _getStringIndexFromMatrix(row,col,ncol2))+=value(matrix1, _getStringIndexFromMatrix(row,colRes,n))*value(matrix2, _getStringIndexFromMatrix(colRes,col,ncol2));
				}
			}
	}
	return result;
}

/*
.Function._matricialSum:
..summary:Computes the matricial sum between two matrixes
..signature:_matricialSum(matrix1,matrix2,ncol)
..param.matrix1:The first matrix.
...type:String<TAlphabet, TSpec1>&
..param.matrix2:The second matrix.
...type:String<TAlphabet, TSpec1>&
..param.ncol:The number of columns of both the matrices.
...type:nolink:unsigned int
..returns:The sum of the two matrices (another matrix).
..remarks:The number of rows and columns of matrix1 must be equal to the number of rows and columns of matrix2.
*/

template <typename TAlphabet, typename TSpec1, typename TSpec2>
String<TAlphabet> _matricialSum(String<TAlphabet, TSpec1>& matrix1,
					  String<TAlphabet, TSpec2>& matrix2, unsigned int ncol)
{
	unsigned int nrow=length(matrix1)/ncol;
	String<TAlphabet> result;
	unsigned int lenghtnew=ncol*nrow;
	fill(result,lenghtnew,0);

	for(unsigned int row = 0; row < nrow; row++){
			for(unsigned int col = 0; col < ncol; col++){
					value(result, _getStringIndexFromMatrix(row,col,ncol))=value(matrix1, _getStringIndexFromMatrix(row,col,ncol))+ value(matrix2, _getStringIndexFromMatrix(row,col,ncol));
			}
	}
	return result;
}

//////////////////////////////////////////////////////////////////////////////
// _matricialDifference
//////////////////////////////////////////////////////////////////////////////

/*
.Function._matricialDifference:
..summary:Computes the matricial difference between two matrixes
..signature:_matricialDifference(matrix1,matrix2,n)
..param.matrix1:The first matrix.
...type:String<TAlphabet, TSpec1>&
..param.matrix2:The second matrix.
...type:String<TAlphabet, TSpec1>&
..param.ncol:The number of columns of both the matrices.
...type:nolink:unsigned int
..returns:The difference of the two matrices (another matrix).
..remarks:The number of rows and columns of matrix1 must be equal to the number of rows and columns of matrix2.
*/

template <typename TAlphabet, typename TSpec1, typename TSpec2>
String<TAlphabet> _matricialDifference(String<TAlphabet, TSpec1>& matrix1,
					  String<TAlphabet, TSpec2>& matrix2, unsigned int ncol)
{
	unsigned int nrow=length(matrix1)/ncol;
	String<TAlphabet> result;
	unsigned int lenghtnew=ncol*nrow;
	resize(result,lenghtnew);

	for(unsigned int row = 0; row < nrow; row++){
			for(unsigned int col = 0; col < ncol; col++){
					value(result, _getStringIndexFromMatrix(row,col,ncol))=value(matrix1, _getStringIndexFromMatrix(row,col,ncol))- value(matrix2, _getStringIndexFromMatrix(row,col,ncol));
			}
	}

	return result;
}

//////////////////////////////////////////////////////////////////////////////
// _computeInverseMatrix
//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeInverseMatrix:
..summary:Computes the inverse matrix of a given matrix
..signature:_computeInverseMatrix(matrix,n)
..param.matrix:The matrix in input.
...type:String<TAlphabet, TSpec1>&
..param.n:The number of columns of the matrix.
...type:nolink:unsigned int
..returns:The inverse matrix of the matrix.
*/

template <typename TAlphabet, typename TSpec1>
String<TAlphabet> _computeInverseMatrix(String<TAlphabet, TSpec1>& matrix, unsigned int n)
{
	String<TAlphabet> result;
	unsigned int lengthnew=n*n;
	fill(result,lengthnew,0);

	//copy the matrix in result, since the procedure is in-place
	String<TAlphabet> tmp = matrix;

	//lu decomposition of a in-place
	String<TAlphabet> indx=_ludcmp(tmp,n);

	String<TAlphabet> col;

	unsigned int i;

	// inverse by columns
	for (unsigned int j=0; j<n; j++) {
		fill(col,n,0);
		if(j>0){
			for(i=0; i<n; i++){
				value(col,i)=0;
			}
		}
		value(col,j) = 1;

		_lubksb(tmp,n,indx,col);

		for (i=0; i<n; i++) {
			value(result, _getStringIndexFromMatrix(i,j,n))= value(col,i);
		}
	}
	return result;
}


template <typename TAlphabet>
void _lubksb(String<TAlphabet> &a, int n, String<TAlphabet> &indx, String<TAlphabet>  &b)
{
  int i, ii=0,ip,j;
  double sum;

  for (i=1; i<=n; i++) {
    ip = value(indx,i-1);
    sum = value(b,ip-1);
    value(b,ip-1) = value(b,i-1);
    if (ii) {
      for (j=ii;j<=i-1;j++)
	  {
		  sum -=value(a,_getStringIndexFromMatrix(i-1,j-1,n))*value(b,j-1);
	  }
    }
    else
		if (sum){
			ii=i;
		}
    value(b,i-1) = sum;
  }
  for (i=n; i>=1; i--) {
    sum = value(b,i-1);
	for (j=i+1; j<=n; j++){
		sum -= value(a,_getStringIndexFromMatrix(i-1,j-1,n))*value(b,j-1);
	}
    value(b,i-1) = sum/value(a,_getStringIndexFromMatrix(i-1,i-1,n));
  }
}

#define TINY 1.0e-20
template <typename TAlphabet>
String<TAlphabet> _ludcmp(String<TAlphabet> &result, int n)
{
  int i, imax, j, k,d;
  double big,dum,sum,temp;
  String<TAlphabet> vv;
  fill(vv, n, 1.0);


  d = 1.0;
  for (i=1; i<=n; i++) {
	big = 0.0;
	for (j=1; j<=n; j++){
		if ((temp=fabs(value(result, _getStringIndexFromMatrix(i-1,j-1,n))))>big){
			big = temp;
		}
	}
	if (big==0.0) {
		std::cout<<"Singular matrix in routine ludcmp" << std::endl;
		exit(1);
	}

	value(vv, _getStringIndexFromMatrix(0,i-1,n)) = 1.0/big;
  }
  String<TAlphabet> indx;
  resize(indx,n);
  for (j=1; j<=n; j++) {
    for (i=1; i<j; i++) {
      sum = value(result, _getStringIndexFromMatrix(i-1,j-1,n));
	  for (k=1; k<i; k++) {
		  sum -= value(result, _getStringIndexFromMatrix(i-1,k-1,n))*value(result, _getStringIndexFromMatrix(k-1,j-1,n));
	  }
	  value(result, _getStringIndexFromMatrix(i-1,j-1,n)) = sum;
    }
    big = 0.0;
    for (i=j; i<=n; i++) {
      sum = value(result, _getStringIndexFromMatrix(i-1,j-1,n));
	  for (k=1; k<j; k++){
		sum -= value(result, _getStringIndexFromMatrix(i-1,k-1,n))*value(result, _getStringIndexFromMatrix(k-1,j-1,n));
	  }
  	  value(result, _getStringIndexFromMatrix(i-1,j-1,n)) = sum;
      if ((dum = value(vv, _getStringIndexFromMatrix(0,i-1,n))*fabs(sum))>=big) {
		big = dum;
		imax = i;
      }
    }
    if (j != imax) {
      for (k=1; k<=n; k++) {
		dum = value(result, _getStringIndexFromMatrix(imax-1,k-1,n));
		value(result, _getStringIndexFromMatrix(imax-1,k-1,n)) = value(result, _getStringIndexFromMatrix(j-1,k-1,n));
		value(result, _getStringIndexFromMatrix(j-1,k-1,n)) = dum;
      }
      d = -(d);
	  value(vv, _getStringIndexFromMatrix(0,imax-1,n))=value(vv, _getStringIndexFromMatrix(0,j-1,n));
    }

	value(indx, _getStringIndexFromMatrix(0,j-1,n)) = imax;

	if (value(result, _getStringIndexFromMatrix(j-1,j-1,n)) == 0.0){
		value(result, _getStringIndexFromMatrix(j-1,j-1,n)) = TINY;
	}
    if (j!=n) {
      dum = 1.0/(value(result, _getStringIndexFromMatrix(j-1,j-1,n)));
	  for (i=j+1; i<=n; i++){
		  value(result, _getStringIndexFromMatrix(i-1,j-1,n)) *= dum;
	  }
    }
  }

 return indx;
}

//AUXILIARY CODE

//////////////////////////////////////////////////////////////////////////////
// printAllLetters
//////////////////////////////////////////////////////////////////////////////

/*
	Only a test function
*/

/*
template<typename TAlphabet>
void
printAllLetters(String<TAlphabet>& str1) {
	typedef typename Iterator<String<TAlphabet> >::Type TIter;
	TIter itStr = begin(str1);
	TIter itStrEnd = end(str1);
	for(;itStr != itStrEnd; ++itStr)
	{
		std::cout << value(itStr) << ',';
	}
	std::cout << std::endl;
}

template<typename TAlphabet>
void
printAllLetters(String<TAlphabet>& str1, int n) {
	typedef typename Iterator<String<TAlphabet> >::Type TIter;
	TIter itStr = begin(str1);
	TIter itStrEnd = end(str1);
	int count=1;
	for(;itStr != itStrEnd; ++itStr)
	{
		std::cout << value(itStr) << ',';
		if ((count%n)==0)
		{
				std::cout << std::endl;
		}
		count++;
	}
	std::cout << std::endl;
}

template<typename TAlphabet>
void
printAllLettersFile(String<TAlphabet>& str1, int n, FILE *f) {
	typedef typename Iterator<String<TAlphabet> >::Type TIter;
	TIter itStr = begin(str1);
	TIter itStrEnd = end(str1);
	int count=1;
	for(;itStr != itStrEnd; ++itStr)
	{
		fprintf(f,"%lf ",value(itStr));
		if ((count%n)==0)
		{
				fprintf(f,"\n");
		}
		count++;
	}
}
*/

//////////////////////////////////////////////////////////////////////////////
// printMarkovModel
//////////////////////////////////////////////////////////////////////////////

/*
	Only a test function
*/

/*
template <typename TAlphabet, typename TFloat, typename TSpec>
void printMarkovModel(MarkovModel<TAlphabet, TFloat, TSpec> & mm)
{
	unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
	unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) mm.order);

	for(unsigned int col = 0; col < column_size;++col) {
		String<TAlphabet> qgram;
		unhash(qgram, col, mm.order);
		std::cout << qgram << ",";
	}
	std::cout << std::endl;
	for(unsigned int row = 0; row < column_size;++row) {
		String<TAlphabet> qgram;
		unhash(qgram, row, mm.order);
		std::cout << qgram << ",";
		for(unsigned int col = 0; col < column_size;++col) {
			std::cout << value(mm.transition, row * column_size + col) << ",";
		}
		std::cout << std::endl;
	}

}
*/

//////////////////////////////////////////////////////////////////////////////
// Interface
//////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TFloat, typename TSpec, typename TStringSet>
void buildMarkovModel(MarkovModel<TAlphabet, TFloat, TSpec> & mm,
					  TStringSet & strings)
{
	mm.build(strings);
}



///////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TFloat, typename TSpec>
void setMarkovModel(MarkovModel<TAlphabet, TFloat, TSpec> & mm,
					String<TFloat> &iTransition)
{
	mm.set(iTransition);
}



template <typename TAlphabet, typename TFloat, typename TSpec>
void setMarkovModel(MarkovModel<TAlphabet, TFloat, TSpec> & mm,
					String<TFloat> &iTransition, 
					String<TFloat> &iStationaryDistribution)
{
	mm.set(iTransition, iStationaryDistribution);
}

///////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TFloat, typename TSpec, typename TStringSet>
TFloat emittedProbability(MarkovModel<TAlphabet, TFloat, TSpec> & mm,
						  TStringSet &string)
{
	return mm.emittedProbability(string);
}


template <typename TAlphabet, typename TFloat, typename TSpec>
TFloat emittedProbability(MarkovModel<TAlphabet, TFloat, TSpec> & mm,
						  String<TAlphabet> &string)
{
	return mm.emittedProbability(string);
}


///////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TFloat, typename TSpec>
void write(FILE *file, 
		   MarkovModel<TAlphabet, TFloat, TSpec> & mm )
{
	mm.write(file);
}


///////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TFloat, typename TSpec>
void read(FILE *file, 
		  MarkovModel<TAlphabet, TFloat, TSpec> & mm )
{
	mm.read(file);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace seqan

#endif //#ifndef SEQAN_HEADER_...
