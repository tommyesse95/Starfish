/* *****************************************************
 * (c) 2012 Particle In Cell Consulting LLC
 * 
 * This document is subject to the license specified in 
 * Starfish.java and the LICENSE file
 * *****************************************************/

package starfish.core.solver;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import starfish.core.common.Starfish;
import starfish.core.common.Starfish.Log;

/*this class provides support for a sparse matrix*/
public class Matrix 
{
    /*variables*/
    protected HashMap<Integer,Double> data[];	//data for each row consisting of column / value
    //protected ArrayList<Double> value[];    /*data held by the matrix*/
    //protected ArrayList<Integer> jx[];	    /*j index in full matrix for each value*/
    public final int nr;		    /*number of rows and values in each row*/
    	
    /**sparse matrix constructor */
    public Matrix(int nr)
    {
	this.nr = nr;	    //save number of rows
	
	data = new HashMap[nr];
	
	for (int i=0;i<nr;i++)
	    data[i] = new HashMap<Integer,Double>();
    }

    /**copy constructor */
    public static Matrix copy(Matrix A)
    {
	Matrix C = new Matrix(A.nr);
	
	for (int i=0;i<A.nr;i++)
	    C.data[i] = (HashMap<Integer,Double>)A.data[i].clone();	    
	return C;
    }

    
    /**clears (sets to zero) a single row*/
    public void clearRow(int i)
    {
	data[i] = new HashMap<Integer,Double>();
    }

    /**returns the value held by full matrix at row i and column j*/
    public double get(int i, int j)
    {
	Double val = data[i].get(j);	//returns null if not found, so need object
	if (val==null) return 0; else return val;
    }

    /**sets value at row i, column j in full matrix*/
    public void set(int i, int j, double val)
    {
	data[i].put(j, val);	
    }

    /**add value to row r, column c in full matrix*/
    public void add(int i, int j, double val)
    {
	set(i,j, get(i,j)+val);
    }

    /** copies single row between matrixes*/
    public void copyRow(Matrix A, int i)
    {
	assert(nr==A.nr);	
	data[i] = (HashMap<Integer,Double>)A.data[i].clone();	
    }
    
    /**add value to row r, column c in full matrix*/
    public void subtract(int i, int j,  double val)
    {
	add(i,j,-val);
    }
    
    /**returns A-B*/
    public Matrix subtract(Matrix B)
    {
	assert(nr==B.nr);
	
	Matrix R = copy(this);
		
	for (int i=0;i<nr;i++)
	{
	    for(Map.Entry<Integer, Double> it : data[i].entrySet())
	    {
		int j = it.getKey();
		double val = it.getValue();
		R.add(i,j,-val);
	    }
	}
	return R;
    }


    /**returns A-diag(B), for now defined only for identical matrices (nv is equal)*/
    public Matrix subtractDiag(double b[])
    {
	assert(nr==b.length);
	
	Matrix R = Matrix.copy(this);
	
	for (int i=0;i<nr;i++)
	    R.set(i, i, get(i,i)-b[i]);
		 
	return R;
    }

    /** performs matrix matrix multiplication
     * 
     * @param A matrix to multiply by
     * @return R=M*A
     */
    public Matrix mult(Matrix A)
    {
	Matrix R = new Matrix(nr);
	
	/*simple algorithm, not optimized in any way*/
	for (int ri=0;ri<nr;ri++)
	    for (int rj=0;rj<nr;rj++)
	    {
		/*multiply rj in M by ci in A*/
		double prod=0;
		for (int k=0;k<nr;k++)
		{
		    prod += get(ri,k)*A.get(k,rj);
		}
		R.set(ri, rj, prod);
	    }
	   	
	return R;
    }

      /**performs matrix vector multiplication
    @return A*x*/
    public double[] mult(double x[])
    {
	double b[] = new double[this.nr];
	mult(x,b);
	return b;
    }
    
    /**performs matrix vector multiplication and stores it in result vector
    @return A*x*/
    public void mult(double x[], double result[])
    {
	for (int r=0;r<nr;r++)
	{
	    double prod=0;
	     for (Map.Entry<Integer, Double> it : data[r].entrySet())
	     {
		 int j = it.getKey();
		 double val = it.getValue();
		 prod += val * x[j];
	     }
	    result[r] = prod;	
	}    	
    }

    /**multiplies value held in full matrix row i, column j, by value*/
    public void mult(int i, int j,  double val)
    {
	data[i].put(j, data[i].get(j)*val);
    }

    /**multiplies entire row by s*/
    public void multRow(int i, double s)
    {
	for (Map.Entry<Integer, Double> it : data[i].entrySet())
	{
	    int j = it.getKey();
	    double val = it.getValue();
	    it.setValue(val*s);
	}
    }

    /**multiplies one row of the matrix by a vector but
     * excludes the diagonal element  */
    public double multRowNonDiag(double x[], int i)
    {
	double prod=0;
	
	for (Map.Entry<Integer, Double> it : data[i].entrySet())
	{
	    int j = it.getKey();
	    double val = it.getValue();
	    if (j!=i) prod+=val*x[j];
	}
	
	return prod;
    }
	
    /**returns the identity matrix of size of A*/
    public Matrix identity()
    {
	Matrix I = new Matrix(nr);		/*diagonal matrix*/
        
	for (int i=0;i<nr;i++)
	    I.set(i, i, 1);

	return I;
    }

    /**returns a new matrix which is the diagonal of the specified one*/
    public Matrix diag_matrix()
    {
	Matrix D = new Matrix(nr);		/*diagonal matrix*/

	for (int i=0;i<D.nr;i++)
	    D.set(i, i, get(i,i));

	return D;
    }

    /**returns a vector containing the diagonal*/
    public double[]diag(Matrix A)
    {
	double D[]= new double[A.nr];		/*diagonal matrix*/

	for (int i=0;i<A.nr;i++)
	    D[i]=A.get(i,i);

	return D;
    }

    /**returns the inverse
     * NOTE: this is currently defined only for diagonal matrix!
     */
    public Matrix inverse()
    {
	//make sure we have a diagonal matrix*/
	for (int i=0;i<nr;i++)
	    if (data[i].size()>1) Log.error("Matrix inverse currently defined only for diagonal matrixes");
	
	Matrix I = new Matrix(nr);		/*diagonal matrix*/

	for (int i=0;i<I.nr;i++)
	    I.set(i,i, 1.0/get(i,i));	/*inverse of a diagonal is 1/D[i,i]*/

	return I;
    }

    /**returns the transpose
     */
    public Matrix transpose()
    {
	Matrix T = new Matrix(nr);		/*diagonal matrix*/

	for (int i=0;i<nr;i++)
	    for (Map.Entry<Integer, Double> it : data[i].entrySet())
	    {
		int j = it.getKey();
		double val = it.getValue();
		T.set(j,i,val);
	    }
	return T;		
    }
	
 

    /**prints self*/
    public void print()
    {
	for (int i=0;i<nr;i++)
	{
	    for (int j=0;j<nr;j++)
	    {
		System.out.printf("%8.2g", get(i,j));
	    }
	    System.out.println();
	}
	System.out.println();
    }	
    
    /**creates 2D transformation matrix
     * M = T*R*S 
     */
    public static Matrix makeTransformationMatrix(double scaling[], double theta, double translation[])
    {
	/*convert to radians*/
	theta *= Math.PI/180.0;
	
	/*rotation matrix
	 * 
	 * R = [C -S  0]
	 *     [S  C  0]
	 *     [0  0  1]
	 */
	Matrix R = new Matrix(3);
	R.set(0, 0, Math.cos(theta));
	R.set(1, 0, -Math.sin(theta));
	R.set(0, 1, Math.sin(theta));
	R.set(1, 1, Math.cos(theta));
	R.set(2, 2, 1);
	
	/*scaling matrix
	 *
	 * S = [sx 0  0]
	 *     [0  sy 0]
	 *     [0  0  1]
	 */
	Matrix S = new Matrix(3);
	S.set(0,0,scaling[0]);
	S.set(1,1,scaling[1]);
	S.set(2,2,1);
	
	/*translation matrix
	 *
	 * T = [1 0 vx]
	 *     [0 1 vy]
	 *     [0 0  0]
	 */
	Matrix T = new Matrix(3);
	T.set(0, 0, 1);
	T.set(0, 2, translation[0]);
	T.set(1, 1, 1);
	T.set(1, 2, translation[1]);
	T.set(2, 2, 1);

	/*apply transformation, rotation->scaling->translation*/
	return T.mult(R.mult(S));
    }
    
    /*cleans up the matrix by removing zero entries*/
    void removeZeros() 
    {
	for (int i = 0;i<nr;i++)
	{
	    //new empty data
	    HashMap<Integer,Double> d = new HashMap<Integer,Double>();
	    
	    for (Map.Entry<Integer, Double> it : data[i].entrySet())
	    {
		int j = it.getKey();
		double val = it.getValue();
		if (val!=0) d.put(j, val);
	    }
	    data[i] = d;
	}
    }
    
    
    /*computes LU decomposition of the matrix without pivoting
    based on the algorithm in Numerical Analysis
    */
    public Matrix[] decomposeLU() throws UnsupportedOperationException
    {
	Log.message("Computing LU decomposition, nr = "+nr);
	Matrix L = new Matrix(nr);
	Matrix U = new Matrix(nr);
	Matrix A = this;    //to make the eq below more explicit
	double s;
	
	//step 1
	L.set(0,0,1);
	U.set(0,0,A.get(0,0));
	if (A.get(0,0)==0) throw new UnsupportedOperationException(); 
	
	//step 2
	for (int j=1;j<nr;j++)
	{
	    U.set(0,j,A.get(0,j)/L.get(0,0));
	    L.set(j,0,A.get(j,0)/U.get(0,0));
	}
	
	//step 3
	for (int i=1;i<nr-1;i++)
	{
	    System.out.printf("%d of %d\n",i,nr);
	    //step 4
	    L.set(i,i,1);
	    s = A.get(i,i);
	    for (int k=0;k<i;k++)
		s-=L.get(i,k)*U.get(k,i);
	    U.set(i,i,s);
	    if (s==0) throw new UnsupportedOperationException();
	    
	    //step 5
	    double Lii = L.get(i,i);
	    double Uii = U.get(i,i);
	    
	    for (int j=i+1;j<nr;j++)
	    {
		s = A.get(i,j);
		for (int k=0;k<i;k++)
		    s -= L.get(i,k)*U.get(k,j);
		U.set(i, j, s/Lii);
		
		s = A.get(j,i);
		for (int k=0;k<i;k++)
		    s -= L.get(j,k)*U.get(k,i);
		L.set(j,i,s/Uii);
	    }
	}
	
	//step 6
	L.set(nr-1,nr-1,1);
	s = A.get(nr-1,nr-1);
	for (int k=0;k<nr-1;k++)
	    s -= L.get(nr-1,k)*U.get(k,nr-1);
	U.set(nr-1,nr-1,s);
	if (s==0) throw new UnsupportedOperationException();
	
	L.removeZeros();
	U.removeZeros();
	
	Log.debug("LU decomposition complete");
	Matrix ret[] = {L,U};
	/*A.print();
	System.out.println();
	L.print();
	System.out.println();
	U.print();*/
	return ret;	
    }
    
    
}
