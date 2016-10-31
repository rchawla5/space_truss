/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-15, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Contains CTruss class implementation.

*********************************************/
#include <iomanip>
#include <sstream>
#include "truss.h"
#include "matlib.h"
#include "matrixtoolbox.h"
#include "..\library\fileio.h"
#include "..\library\parser.h"


const int MAXCHARS = 80;
std::string szInputString;
std::string szComment ("**");

/* ==================================================================
   ======================= CTruss class =============================
   ================================================================== */

CTruss::CTruss ()
// ---------------------------------------------------------------------------
// Function: default ctor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nNodes = m_nElements = m_nDOF = 0;
    m_nDebugLevel = m_nLineNumber = 0;
}

CTruss::~CTruss ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CTruss::Banner (std::ostream& OF) const
// ---------------------------------------------------------------------------
// Function: prints the program banner on the output stream
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    OF << '\n';
    OF << "\t\t--------------------------------------------" << '\n';
    OF << "\t\t        Planar Truss Analysis Program       " << '\n';
    OF << "\t\tIntroduction to Structural Analysis & Design" << '\n';
    OF << "\t\t          (c) 2000-15, S. D. Rajan          " << '\n';
    OF << "\t\t        Enhanced By: Rashmi Chawla         " << '\n';
    OF << "\t\t--------------------------------------------" << '\n';
}

void CTruss::PrepareIO (int argc, char* argv[])
// ---------------------------------------------------------------------------
// Function: opens the input and output files
// Input:    command line arguments (currently unused)
// Output:   none
// ---------------------------------------------------------------------------
{
    if (argc == 1)
    {
        // open the input file
        OpenInputFileByName ("Complete input file name: ", m_FileInput,
                             std::ios::in);

        // open the output file
        OpenOutputFileByName ("Complete output file name: ", m_FileOutput,
                              std::ios::out);
    }
    else if (argc == 3) // spacetruss input_file output_file
    {
        m_FileInput.open (argv[1], std::ios::in);
        if (!m_FileInput)
            ErrorHandler (ERRORCODE::CANNOTOPENIFILE);
        m_FileOutput.open (argv[2], std::ios::out);
        if (!m_FileOutput)
            ErrorHandler (ERRORCODE::CANNOTOPENOFILE);
        std::cout << "\n";
        std::cout << argv[1] << " opened as input file.\n";
        std::cout << argv[2] << " opened as output file.\n";
    }
	else
    {
        ErrorHandler (ERRORCODE::INVALIDCOMMANDLINE);
    }

    // print banner
    Banner (m_FileOutput);
}

void CTruss::ReadProblemSize ()
// ---------------------------------------------------------------------------
// Function: Reads the size of the problem being solved
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    CParser Parse;

    // header line
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        ErrorHandler (ERRORCODE::INVALIDINPUT);
    // read the problem description
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        ErrorHandler (ERRORCODE::INVALIDINPUT);
    // header line
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        ErrorHandler (ERRORCODE::INVALIDINPUT);
    // read number of nodes, elements and debug level
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        ErrorHandler (ERRORCODE::INVALIDINPUT);
    std::istringstream szFormatString (szInputString);
    szFormatString >> m_nNodes >> m_nElements >> m_nDebugLevel;
    if (szFormatString.fail() || szFormatString.bad())
        ErrorHandler (ERRORCODE::INVALIDINPUT);

    // check data for validity
    if (m_nNodes <= 1) ErrorHandler (ERRORCODE::INVALIDNUMNODES);
    if (m_nElements <= 0) ErrorHandler (ERRORCODE::INVALIDNUMELEMENTS);
    if (m_nDebugLevel < 0 || m_nDebugLevel > 1) ErrorHandler (ERRORCODE::INVALIDDEBUGCODE);

    // dynamic memory allocations for arrays
    SetSize ();
}

void CTruss::SetSize ()
// ---------------------------------------------------------------------------
// Function: Carries out memory allocation for all the major arrays
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // allocate space for nodal data
    m_NodalData.SetSize (m_nNodes);
    // allocate space for nodal response data
    m_NodalResponseData.SetSize (m_nNodes);
    // allocate space for element data
    m_ElementData.SetSize (m_nElements);
    // allocate space for element response data
    m_ElementResponseData.SetSize (m_nElements);

    // allocate and initialize major matrices
    m_nDOF = 3*m_nNodes;
    m_dSSM.SetSize (m_nDOF, m_nDOF);
    m_dSND.SetSize (m_nDOF, 1);
    m_dSNF.SetSize (m_nDOF, 1);
    // initialize all elements of matrices
	m_dSSM.Set (0.0);
    m_dSND.Set (0.0);
    m_dSNF.Set (0.0);
}

void CTruss::ReadTrussModel ()
// ---------------------------------------------------------------------------
// Function: Read the truss model data from the input file
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // read problem size
    ReadProblemSize ();
    
    int i, nN;
    float fX, fY,fZ;
    CParser Parse;
    
    
	// header line
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        ErrorHandler (ERRORCODE::INVALIDINPUT);
    
	
	// read nodal coordinates
    for (i=1; i <= m_nNodes; i++)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            ErrorHandler (ERRORCODE::INVALIDINPUT);
        else
        {
            std::istringstream szFormatString (szInputString);
            szFormatString >> nN >> fX >> fY  >> fZ;
            if (szFormatString.fail() || szFormatString.bad())
                ErrorHandler (ERRORCODE::INVALIDINPUT);
        }
        if (nN <= 0 || nN > m_nNodes) 
            ErrorHandler (ERRORCODE::INVALIDNODENUM);
        m_NodalData(nN).SetCoords (fX, fY, fZ);
    }

    
	
	// header line
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        ErrorHandler (ERRORCODE::INVALIDINPUT);
    
	
	
	// read nodal fixity conditions
    int nXFC, nYFC, nZFC;
    for (i=1; i <= m_nNodes; i++)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
            ErrorHandler (ERRORCODE::INVALIDINPUT);
        else
        {
            std::istringstream szFormatString (szInputString);
            szFormatString >> nN >> nXFC >> nYFC >> nZFC;
            if (szFormatString.fail() || szFormatString.bad())
                ErrorHandler (ERRORCODE::INVALIDINPUT);
        }
        if (nN <= 0 || nN > m_nNodes) ErrorHandler (ERRORCODE::INVALIDNODENUM);
        if (nXFC < 0 || nXFC > 1) ErrorHandler (ERRORCODE::INVALIDNODALFIXITY);
        if (nYFC < 0 || nYFC > 1) ErrorHandler (ERRORCODE::INVALIDNODALFIXITY);
		if (nZFC < 0 || nZFC > 1) ErrorHandler (ERRORCODE::INVALIDNODALFIXITY);
        m_NodalData(nN).SetFixity (nXFC, nYFC, nZFC);
    }
	
	
	// header line
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        ErrorHandler (ERRORCODE::INVALIDINPUT);
    
	
	// read nodal loads
    float fXF, fYF, fZF;
    for (i=1; i <= m_nNodes; i++)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
            ErrorHandler (ERRORCODE::INVALIDINPUT);
        else
        {
            std::istringstream szFormatString (szInputString);
            szFormatString >> nN >> fXF >> fYF >> fZF;
            if (szFormatString.fail() || szFormatString.bad())
                ErrorHandler (ERRORCODE::INVALIDINPUT);
        }
        if (nN <= 0 || nN > m_nNodes) ErrorHandler (ERRORCODE::INVALIDNODENUM);
        //
        m_NodalData(nN).SetLoads (fXF, fYF, fZF);
    }
	
	
   
	// header line
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        ErrorHandler (ERRORCODE::INVALIDINPUT);
    
	
	// read element data
    int nE, nSN, nEN;
    float fA, fE;
    for (i=1; i <= m_nElements; i++)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
            ErrorHandler (ERRORCODE::INVALIDINPUT);
        else
        {
            std::istringstream szFormatString (szInputString);
            szFormatString >> nE >> nSN >> nEN >> fA >> fE;
            if (szFormatString.fail() || szFormatString.bad())
                ErrorHandler (ERRORCODE::INVALIDINPUT);
        }
        if (nE <= 0 || nE > m_nElements) ErrorHandler (ERRORCODE::INVALIDELEMENTNUM);
        if (nSN <= 0 || nSN > m_nNodes) ErrorHandler (ERRORCODE::INVALIDNODENUM);
        if (nEN <= 0 || nEN > m_nNodes) ErrorHandler (ERRORCODE::INVALIDNODENUM);
        if (fA <= 0.0f) ErrorHandler (ERRORCODE::INVALIDCSAREA);
        if (fE <= 0.0f) ErrorHandler (ERRORCODE::INVALIDYM);
        m_ElementData(i).SetData (nSN, nEN, fA, fE);
    }
	
	
	// header line
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        ErrorHandler (ERRORCODE::INVALIDINPUT);
	
	
	// read elemental temperature data
	int nET;
	float falpha, fCT;
	for (i=1; i <= m_nElements; i++)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
            ErrorHandler (ERRORCODE::INVALIDINPUT);
        else
        {
            std::istringstream szFormatString (szInputString);
            szFormatString >> nET >> falpha >> fCT;
            if (szFormatString.fail() || szFormatString.bad())
                ErrorHandler (ERRORCODE::INVALIDINPUT);
        }
        //if (nE <= 0 || nE > m_nElements) ErrorHandler (ERRORCODE::INVALIDELEMENTNUM);
        //if (nSN <= 0 || nSN > m_nNodes) ErrorHandler (ERRORCODE::INVALIDNODENUM);
        //if (nEN <= 0 || nEN > m_nNodes) ErrorHandler (ERRORCODE::INVALIDNODENUM);
        //if (fA <= 0.0f) ErrorHandler (ERRORCODE::INVALIDCSAREA);
        //if (fE <= 0.0f) ErrorHandler (ERRORCODE::INVALIDYM);
        m_ElementData(i).SetTemperatureData(falpha,fCT);
    }

     // end keyword?
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        ErrorHandler (ERRORCODE::INVALIDINPUT);
    if (szInputString.substr(0,4) != "*end")
        ErrorHandler (ERRORCODE::MISSINGEND);

    // construct structural nodal load vector due to external load
	for (i=1; i <= m_nNodes; i++)
    {
        m_NodalData(i).GetLoads (fXF, fYF, fZF);
        m_dSNF(3*i-2,1) = fXF;
        m_dSNF(3*i-1,1) = fYF;
		m_dSNF(3*i,1) = fZF;
    }
}

void CTruss::Analyze ()
// ---------------------------------------------------------------------------
// Function: Implements the FEA steps
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
	// construct structural stiffness matrix
	ConstructK ();

	// impose boundary conditions
	ImposeBC ();

	// solve for the nodal displacements
	Solve ();

	// compute element response
	Response ();
	
    // create output file
	CreateOutput ();
}

void CTruss::ConstructK ()
// ---------------------------------------------------------------------------
// Function: Constructs the structural stiffness matrix 
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    int i, j, k;
    CVector<int> nE(6);
    CMatrix<double> dkl(2,2), dkg(6,6);
    CMatrix<double> dT(2,6), dTT(6,2);
    CMatrix<double> dTemp(2,6);
    float fX1, fX2, fY1, fY2, fZ1, fZ2, fL;
    float fA, fE;
    int nSN, nEN;
    
    // initialize
    dT.Set (0.0);
    
    // loop thro' all elements
    for (i=1; i <= m_nElements; i++)
    {
        // construct k(4x4) in two steps
        m_ElementData(i).GetData (nSN, nEN, fA, fE);
        m_NodalData(nSN).GetCoords (fX1, fY1, fZ1);
        m_NodalData(nEN).GetCoords (fX2, fY2, fZ2);
        fL = sqrt((fX2-fX1)*(fX2-fX1) +
                  (fY2-fY1)*(fY2-fY1) +
				  (fZ2-fZ1)*(fZ2-fZ1));
        // form AE/L 
        double dAEOL = double(fA*fE/fL);
        // construct klocal
        dkl(1,1) = dAEOL;	dkl(1,2) = -dAEOL;	
        dkl(2,1) = -dAEOL;	dkl(2,2) = dAEOL;	
		
        // local-to-global transformation matrix
        double dl = double((fX2-fX1)/fL);
        double dm = double((fY2-fY1)/fL);
		double dn = double((fZ2-fZ1)/fL);
        dT(1,1) = dl; dT(1,2) = dm; dT(1,3) = dn;
        dT(2,4) = dl; dT(2,5) = dm; dT(2,6) = dn;
        CMatToolBox<double> mattool;
		// transpose of the T matrix
		mattool.Transpose (dT, dTT);
        // form k'*T
        mattool.Multiply (dkl, dT, dTemp);
        // construct kglobal=T(T)*k'*T
        mattool.Multiply (dTT, dTemp, dkg);

        // assemble into structural K
        nE(1) = 3*nSN-2; nE(2) = 3*nSN-1; nE(3) = 3*nSN;
        nE(4) = 3*nEN-2; nE(5) = 3*nEN-1; nE(6) = 3*nEN;
       
		for (j=1; j <= 6; j++)
        {
            int nRow = nE(j);
            for (k=1; k <= 6; k++)
            {
                int nCol = nE(k);
                m_dSSM(nRow, nCol) += dkg(j,k);
            }
        }
			
        // debug?
        if (m_nDebugLevel == 1)
        {
            std::ostringstream szPrompt;
            szPrompt << "Stiffness Matrix for element " << i;
            PrintMatrixRowWise (dkg, szPrompt.str(), m_FileOutput);
        }
    }

    // debug?
    if (m_nDebugLevel == 1)
    {
        PrintMatrixRowWise (m_dSSM, "Structural Stiffness (Before BCs)",
                            m_FileOutput);
    }
}

void CTruss::ImposeBC ()
// ---------------------------------------------------------------------------
// Function: Imposes the essential boundary conditions and temperature stresses
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
// construct structural nodal load due to temperature change
	CVector <int> nL(6);
	float fA, fE;
    int nSN, nEN;
	float falpha,fCT;
	float fX1, fX2, fY1, fY2,fZ1,fZ2, fL;
	CVector <double> q(2),v_dENF(6);
	CMatrix <double> dT(2,6),dTT(6,2);
	q.Set(0.0);
	v_dENF.Set(0.0);
	dT.Set(0.0);
	dTT.Set(0.0);
	for (int i=1; i <= m_nElements; i++)
    {
        m_ElementData(i).GetData (nSN, nEN, fA , fE);
        m_NodalData(nSN).GetCoords (fX1, fY1, fZ1);
        m_NodalData(nEN).GetCoords (fX2, fY2, fZ2);
        fL = sqrt((fX2-fX1)*(fX2-fX1) +
                  (fY2-fY1)*(fY2-fY1) +
				  (fZ2-fZ1)*(fZ2-fZ1));
        // local-to-global transformation matrix
        double dl = double((fX2-fX1)/fL);
        double dm = double((fY2-fY1)/fL);
		double dn = double((fZ2-fZ1)/fL);
        
		dT(1,1) = dl; dT(1,2) = dm, dT(1,3) = dn;
        dT(2,4) = dl; dT(2,5) = dm; dT(2,6) = dn;
		
		m_ElementData(i).GetTemperatureData (falpha,fCT);
		
		q(1) = -fA*fE*falpha*fCT;
		q(2) = fA*fE*falpha*fCT;
		
		CMatToolBox<double> mattool;
		mattool.Transpose(dT,dTT);
		mattool.MatMultVec(dTT,q,v_dENF);
		
		nL(1) = 3*nSN-2; nL(2) = 3*nSN-1; nL(3) = 3*nSN;
        nL(4) = 3*nEN-2; nL(5) = 3*nEN-1; nL(6) = 3*nEN;

		// modification of structural nodal load due to temperature change	
		for (int j=1; j <= 6; j++)
        { 
			m_dSNF(nL(j), 1) = m_dSNF(nL(j), 1) + v_dENF(j);
        }
    }
	
    int nXFC, nYFC, nZFC;
    
    // loop thro' all nodes
    for (int i=1; i <= m_nNodes; i++)
    {
        m_NodalData(i).GetFixity (nXFC, nYFC, nZFC);
        if (nXFC == 1)
        {
            int nGDOF = 3*i-2;
            SuppressDOF (nGDOF);
        }
        if (nYFC == 1)
        {
            int nGDOF = 3*i-1;
            SuppressDOF (nGDOF);
        }
		if (nZFC == 1)
        {
            int nGDOF = 3*i;
            SuppressDOF (nGDOF);
        }
    }
	
	// debug?
    if (m_nDebugLevel == 1)
    {
        PrintMatrixRowWise (m_dSSM, "Structural Stiffness (After BCs)",
                            m_FileOutput);
        PrintMatrixColumnWise (m_dSNF, "Structural Nodal Forces (After BCs)",
                               m_FileOutput);
    }
	
}

 void CTruss::SuppressDOF (const int nEqn)
// ---------------------------------------------------------------------------
// Function: Imposes the essential boundary conditions
// Input:    Equation number
// Output:   none
// ---------------------------------------------------------------------------
{

    for (int j=1; j <= m_nDOF; j++)
    {
        // zero out the row
        m_dSSM(nEqn, j) = 0.0;
        // zero out the column
        m_dSSM(j, nEqn) = 0.0;
    }
    // set diagonal to 1
    m_dSSM(nEqn, nEqn) = 1.0;

    // set RHS to zero
    m_dSNF(nEqn, 1) = 0.0;
}

 void CTruss::Solve ()
// ---------------------------------------------------------------------------
// Function: Solves the system equations for the nodal displacements
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    
    const double TOL = 1.0e-6;
	CVector<double> v_dSND(m_nDOF),v_dSNF(m_nDOF);
	v_dSND.Set(0.0);
	v_dSNF.Set(0.0);
	for (int i = 1; i <= m_nDOF; i++)
	{
		v_dSNF(i) = m_dSNF(i,1);
		v_dSND(i) = m_dSND(i,1);
	}
	CMatToolBox<double> mattool;
    mattool.LDLTFactorization (m_dSSM, TOL);
	mattool.LDLTSolve (m_dSSM,v_dSND,v_dSNF);
	//GaussElimination ( m_dSSM,m_dSND, m_dSNF,TOL);
                   
	for (int i = 1; i <= m_nDOF; i++)
	{
		std::cout << v_dSND(i);
	}	
	for (int i = 1; i <= m_nDOF; i++)
	{
		m_dSND(i,1) = v_dSND(i); 
		m_dSNF(i,1) = v_dSNF(i);
	}
	
	for (int i=1; i <= m_nNodes; i++)
    {
			float m_fXDisp = static_cast<float>(m_dSND(3*i-2,1));
            float m_fYDisp = static_cast<float>(m_dSND(3*i-1,1));
			float m_fZDisp = static_cast<float>(m_dSND(3*i,1));
            m_NodalResponseData(i).SetDisplacements (m_fXDisp, m_fYDisp,m_fZDisp);
   }
 }
 
void CTruss::Response ()
// ---------------------------------------------------------------------------
// Function: Computes the element response
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    int i, j;
    CVector<int> nE(6); // stores the element dof
    CMatrix<double> dT(2,6), dTT(6,2); // transformation matrices
    CMatrix<double> dLD(2,1), dND(6,1); // nodal displacements
    float fX1, fX2, fY1, fY2,fZ1,fZ2, fL;
    float fA, fE;
    int nSN, nEN;
    
    // initialize
    dT.Set (0.0);
    
    // loop thro' all elements
    for (i=1; i <= m_nElements; i++)
    {
        // form strain, stress and force in two steps
        m_ElementData(i).GetData (nSN, nEN, fA, fE);
        m_NodalData(nSN).GetCoords (fX1, fY1, fZ1);
        m_NodalData(nEN).GetCoords (fX2, fY2, fZ2);
        fL = sqrt((fX2-fX1)*(fX2-fX1) +
                  (fY2-fY1)*(fY2-fY1) +
				  (fZ2-fZ1)*(fZ2-fZ1));
        // local-to-global transformation matrix
        double dl = double((fX2-fX1)/fL);
        double dm = double((fY2-fY1)/fL);
		double dn = double((fZ2-fZ1)/fL);
        dT(1,1) = dl; dT(1,2) = dm, dT(1,3) = dn;
        dT(2,4) = dl; dT(2,5) = dm; dT(2,6) = dn;

        // get element nodal displacements
        nE(1) = 3*nSN-2; nE(2) = 3*nSN-1; nE(3) = 3*nSN;
        nE(4) = 3*nEN-2; nE(5) = 3*nEN-1; nE(6) = 3*nEN;
        for (j=1; j <= 6; j++)
        {
            dND(j,1) = m_dSND(nE(j),1);
        }

        // form d'=T*d
        CMatToolBox<double> Mattool;
		Mattool.Multiply ( dT, dND, dLD);

        // strain
        float fStrain = float(dLD(2,1) - dLD(1,1))/fL;
        // stress
        float fStress = fE*fStrain;
        // force
        float fForce = fStress*fA;

        // update with the computed values
        m_ElementResponseData(i).SetData (fStrain,
                                          fStress, fForce);
    }
}

void CTruss::CreateOutput ()
// ---------------------------------------------------------------------------
// Function: Creates the output file containing the results.
//           Currently incomplete.
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    int i;
    // print the problem size
    m_FileOutput << '\n';
    m_FileOutput << "PROBLEM SIZE" << '\n';
    m_FileOutput << "------------" << '\n';
    m_FileOutput << "   Number of nodes : " << m_nNodes << '\n';
    m_FileOutput << "Number of elements : " << m_nElements << '\n';

    // print the nodal coordinates
    m_FileOutput << '\n';
    m_FileOutput << "NODAL COORDINATES" << '\n';
    m_FileOutput << "-----------------" << '\n';
    m_FileOutput << "Node" << "   " << "X-Coordinate"
                           << "    " << "Y-Coordinate" <<"	" << "Z-Coordinate"<< '\n';
    m_FileOutput << "----" << "   " << "------------"
                 << "    " << "------------" << "	" << "-------------" << '\n';
	float fX,fY,fZ;
	for (int i =1 ; i<= m_nNodes; i++)
	{
		m_NodalData(i).GetCoords (fX,fY,fZ);
		m_FileOutput << i <<"		"<< fX <<"		"<< fY <<"		"<< fZ <<"\n";
	}

    // print the nodal fixities
    m_FileOutput << '\n';
    m_FileOutput << "NODAL FIXITIES" << '\n';
    m_FileOutput << "--------------" << '\n';
    m_FileOutput << "Node" << "   " << "X-Fixity"
                           << "    " << "Y-Fixity" << "	"<< "Z-Fixity"<< '\n';
    m_FileOutput << "----" << "   " << "--------"
                 << "    " << "--------" <<"	"<<"---------" << '\n';
	int fFX,fFY,fFZ;
	for (int i =1 ; i<= m_nNodes; i++)
	{
		m_NodalData(i).GetFixity(fFX,fFY,fFZ);
		m_FileOutput << i <<"	"<< fFX <<"		"<< fFY <<"		"<< fFZ <<"\n";
	}

    // print the nodal forces
    m_FileOutput << '\n';
    m_FileOutput << "NODAL FORCES" << '\n';
    m_FileOutput << "------------" << '\n';
    m_FileOutput << "Node" << "   " << "X-Force" 
                           << "    " << "Y-Force" 
						   << "    " << "Z-Force" << '\n';
    m_FileOutput << "----" << "   " << "-------" 
                 << "    " << "-------" << '\n';
	float fXF,fYF,fZF;
	for (int i =1 ; i<= m_nNodes; i++)
	{
		m_NodalData(i).GetLoads (fXF,fYF,fZF);
		m_FileOutput << i <<"	"<< fXF <<"		"<< fYF << "	"<< fZF <<"\n";
	}

    // print the nodal displacements
    m_FileOutput << '\n';
    float fXDisp, fYDisp, fZDisp;
    m_FileOutput << "NODAL DISPLACEMENTS" << '\n';
    m_FileOutput << "-------------------" << '\n';
    m_FileOutput << "Node" << "   " << "X-Displacement"
						   << "    " << "Y-Displacement" 
						   << "		" << "Z-Displacement"
						   <<'\n';
    m_FileOutput << "----" << "   " << "--------------"
						   << "    " << "--------------" 
						   << "    " << "--------------"
						   << '\n';
    for (i=1; i <= m_nNodes; i++)
    {
        m_NodalResponseData(i).GetDisplacements (fXDisp, fYDisp, fZDisp);
        m_FileOutput << std::setw(4) << i << "   "
                     << std::setw(14) << fXDisp 
                     << "    " << std::setw(14) << fYDisp 
					 <<"     " << std::setw(14) << fZDisp << '\n';
    }
    
    // print the element strain, stress and force
    m_FileOutput << '\n';
    m_FileOutput << "ELEMENT RESPONSE" << '\n';
    m_FileOutput << "----------------" << '\n';
    m_FileOutput << std::setw (15)<<"Element" <<"	"<< std::setw (15) << "Strain"<<"	" 
                 << std::setw (15) << "Stress" <<"	"<< std::setw (15) << "Force" << '\n';
    m_FileOutput << "---------------" <<"	"<< "---------------" << "	" 
                 << "---------------" << "  " << "---------------" << '\n';
	for (int i=1; i <= m_nElements; i++)
    {
        float fStrain; // strain
        float fStress; // stress
        float fForce;	 // force
		m_ElementResponseData(i).GetData (fStrain,fStress,fForce);
        m_FileOutput << std::setw (15)<< i <<"	"<< std::setw (15)<< fStrain <<"	"
					 << std::setw (15)<< fStress <<"	"<< std::setw (15)<< fForce<<"\n";
    }
}

void CTruss::TerminateProgram ()
// ---------------------------------------------------------------------------
// Function: Closes input and output files
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{

	// close the input and output files
    m_FileInput.close ();
    m_FileOutput.close ();

    std::cout << "\nExecution completed successfully." << std::endl;
}

void CTruss::ErrorHandler (ERRORCODE nCode) const
// ---------------------------------------------------------------------------
// Function: Displays the error message on the error stream
// Input:    Error #
// Output:   none
// ---------------------------------------------------------------------------
{
    std::cerr << '\n';

    if (nCode == ERRORCODE::INVALIDNUMNODES)                                 // invalid number of nodes
        std::cerr << "Number of nodes must be >= 2.";
    else if (nCode == ERRORCODE::INVALIDNUMELEMENTS)                         // invalid number of elements
        std::cerr << "Number of elements must be >= 1.";
    else if (nCode == ERRORCODE::INVALIDDEBUGCODE)                           // invalid debug level
        std::cerr << "Debug level must be 0 or 1.";
    else if (nCode == ERRORCODE::INVALIDNODENUM)                             // invalid node number
        std::cerr << "Invalid node number";
    else if (nCode == ERRORCODE::INVALIDELEMENTNUM)                          // invalid element number
        std::cerr << "Invalid element number";
    else if (nCode == ERRORCODE::INVALIDCSAREA)                              // invalid x/s area
        std::cerr << "Area must be positive.";
    else if (nCode == ERRORCODE::INVALIDYM)                                  // invalid E
        std::cerr << "Modulus of elasticity must be positive.";
    else if (nCode == ERRORCODE::UNSTABLETRUSS)                              // unstable structure
        std::cerr << "Unstable truss.";
    else if (nCode == ERRORCODE::INVALIDINPUT)                               // invalid input
        std::cerr << "Input file need *end as last line.";
    else if (nCode == ERRORCODE::INVALIDNODALFIXITY)                         // invalid fixity code
        std::cerr << "Nodal fixity code must be 0 or 1.";
    else if (nCode == ERRORCODE::MISSINGEND)                                 // missing end statement
        std::cerr << "Missing *END statement in input file.";
    else if (nCode == ERRORCODE::CANNOTOPENIFILE)                            // cannot open input file
        std::cerr << "Cannot open specified input file.";
    else if (nCode == ERRORCODE::CANNOTOPENOFILE)                            // cannot open output file
        std::cerr << "Cannot open specified output file.";
    else if (nCode == ERRORCODE::INVALIDCOMMANDLINE)                         // need 1 or 3 command line arguments
        std::cerr << "Invalid number of command line arguments.";

    if (nCode != ERRORCODE::UNSTABLETRUSS && nCode != ERRORCODE::INVALIDCOMMANDLINE)
        std::cerr << '\n' << "Error in input file line : " 
                  << m_nLineNumber;

    std::cerr << std::endl;

    // fatal error, exit the program
    exit (1);
}