// $Id: fermact_qprop.cc,v 1.5 2003-12-02 15:44:01 edwards Exp $
/*! \file
 *  \brief Propagator solver for a generic non-preconditioned fermion operator
 *
 *  Solve for the propagator of a generic non-preconditioned fermion operator
 */

#include "chromabase.h"
#include "fermact.h"
#include "actions/ferm/invert/invcg2.h"

using namespace QDP;

//! Propagator of a generic non-preconditioned fermion linear operator
/*! \ingroup qprop
 *
 * This routine is actually generic to all non-preconditioned (not red/black) fermions
 *
 * Compute the lattice fermion for a generic non-red/black fermion
 * using the source in "chi" - so, the source can
 * be of any desired form. The result will appear in "psi", which on input
 * contains an initial guess for the solution.

 * \param psi      quark propagator ( Modify )
 * \param state    gauge field ( Read )
 * \param chi      source ( Read )
 * \param invType  inverter type ( Read (
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

template<typename T>
static 
void qprop_t(const FermionAction<T>& me,
	     T& psi, 
	     const ConnectState& state, 
	     const T& chi, 
	     enum InvType invType,
	     const Real& RsdCG, 
	     int MaxCG, int& ncg_had)
{
  START_CODE("FermionAction::qprop");

  int n_count;
  
  /* Construct the linear operator */
  /* This allocates field for the appropriate action */
  const LinearOperatorProxy<T> A(me.linOp(state));

  switch(invType)
  {
  case CG_INVERTER: 
  {
    /* tmp = M_dag(u) * chi_1 */
    T  tmp;
    A(tmp, chi, MINUS);
    
    /* psi = (M^dag * M)^(-1) chi */
    InvCG2 (A, tmp, psi, RsdCG, MaxCG, n_count);
  }
  break;
  
#if 0
  case MR_INVERTER:
    /* psi = M^(-1) chi */
    InvMR (A, chi, psi, MRover, RsdCG, MaxCG, n_count);
    break;

  case BICG_INVERTER:
    /* psi = M^(-1) chi */
    InvBiCG (A, chi, psi, RsdCG, MaxCG, n_count);
    break;
#endif
  
  default:
    QDP_error_exit("Unknown inverter type", invType);
  }
  
  if ( n_count == MaxCG )
    QDP_error_exit("no convergence in the inverter", n_count);
  
  ncg_had = n_count;
  
  END_CODE("FermionAction::qprop");
}


template<>
void FermionAction<LatticeFermion>::qpropT(LatticeFermion& psi, 
					   const ConnectState& state, 
					   const LatticeFermion& chi, 
					   enum InvType invType,
					   const Real& RsdCG, 
					   int MaxCG, int& ncg_had) const
{
  qprop_t(*this, psi, state, chi, invType, RsdCG, MaxCG, ncg_had);
}

template<>
void FermionAction<LatticeFermion>::qprop(LatticeFermion& psi, 
					  const ConnectState& state, 
					  const LatticeFermion& chi, 
					  enum InvType invType,
					  const Real& RsdCG, 
					  int MaxCG, int& ncg_had) const
{
  qprop_t(*this, psi, state, chi, invType, RsdCG, MaxCG, ncg_had);
}



template<>
void FermionAction<LatticeDWFermion>::qpropT(LatticeDWFermion& psi, 
					     const ConnectState& state, 
					     const LatticeDWFermion& chi, 
					     enum InvType invType,
					     const Real& RsdCG, 
					     int MaxCG, int& ncg_had) const
{
  qprop_t(*this, psi, state, chi, invType, RsdCG, MaxCG, ncg_had);
}

template<>
void FermionAction<LatticeDWFermion>::qprop(LatticeFermion& psi, 
					    const ConnectState& state, 
					    const LatticeFermion& chi, 
					    enum InvType invType,
					    const Real& RsdCG, 
					    int MaxCG, int& ncg_had) const
{
  QDPIO::cerr << "FermionAction<DWF>::qprop - this implementation is empty" << endl;
  QDP_abort(1);
}


