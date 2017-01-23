// -----------------------------------------------------------------
// Calculation of the modenumber with the Chebyshev approximation and the Giusti Luescher
// method. Basic algorithm taken from my implementation in the adjoint QCD code.
// Includes some basic checks with trivial operators to verify the correct approximation of the
// step function.
#include "susy_includes.h"
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Z2 random Twist_Fermion
void rand_TFsourceZ2(Twist_Fermion *src) {
	register int i, j, mu;
	register site *s;
	complex grn;

	// Begin with pure z2 random numbers
	FORALLSITES(i, s)
	{
		clear_TF(&(src[i]));
		for (j = 0; j < DIMF; j++) {                // Site fermions
#ifdef SITERAND
			grn.real = z2_rand_no(&(s->site_prn));
			grn.imag = z2_rand_no(&(s->site_prn));
#else
			grn.real = z2_rand_no(&node_prn);
			grn.imag = z2_rand_no(&node_prn);
#endif
			c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[i].Fsite));
			FORALLDIR(mu)
			{                           // Link fermions
#ifdef SITERAND
				grn.real = z2_rand_no(&(s->site_prn));
				grn.imag = z2_rand_no(&(s->site_prn));
#else
				grn.real = z2_rand_no(&node_prn);
				grn.imag = z2_rand_no(&node_prn);
#endif
				c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[i].Flink[mu]));
			}
			for (mu = 0; mu < NPLAQ; mu++) {         // Plaquette fermions
#ifdef SITERAND
				grn.real = z2_rand_no(&(s->site_prn));
				grn.imag = z2_rand_no(&(s->site_prn));
#else
				grn.real = z2_rand_no(&node_prn);
				grn.imag = z2_rand_no(&node_prn);
#endif
				c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[i].Fplaq[mu]));
			}
		}
	}
}

//#define DEBUG_CG
//#define DEBUG_TEST_CLENSH
//#define DEBUG_TEST_PROJ

// Only for testing:
void trivialOperator(const Real sig, Twist_Fermion *src,
		Twist_Fermion *dest){
	register site* s;
	int i;
	FORALLSITES(i, s)
	{
		scalar_mult_TF(&(src[i]), sig, &(dest[i]));
	}
}

// Simple CG version including shifts
void invertShiftedCG(const Real shift, const Real residgoal,
		const unsigned long maxit, Twist_Fermion *src, Twist_Fermion *dest,
		Real* resid, int* iterations, Twist_Fermion ** workv) {
	*iterations = 0;
	Real normsrc = 1.0;
	*resid = 1.0;
	Twist_Fermion * tmp = workv[0];
	Twist_Fermion* dvec = workv[1];
	Twist_Fermion* rvec = workv[2];
	register site* s;
	int i, ite;
	Real rsq, rsqstop, sumtmp, oldrsq, dvecz, dvecdvec, aval, beta;
	complex ctmp;

#ifdef DEBUG_TEST_PROJ
	trivialOperator(0.2,dest,tmp);
#else
	DSq(dest, tmp);
#endif
	normsrc = 0.0;
	sumtmp = 0.0;
	FORALLSITES(i, s)
	{
		normsrc += (double) magsq_TF(&(src[i]));
		scalar_mult_mult_add_TF(&(tmp[i]), -1.0, &(dest[i]), -shift,
				&(rvec[i]));
		scalar_mult_sum_TF(&(src[i]), 1.0, &(rvec[i])); // dvec=rvec=src-Dsq*dest-shift*dest
		copy_TF(&(rvec[i]), &(dvec[i]));
		sumtmp += (double) magsq_TF(&(rvec[i]));
	}
	g_doublesum(&normsrc);
	if (normsrc < 0.01) { // avoid problems for too small initial vectors
		normsrc = 1.0;
	}
	g_doublesum(&sumtmp);
	rsq = sumtmp;
	rsqstop = residgoal * normsrc;

#ifdef DEBUG_CG
	node0_printf("Initialized Mode CG res=%g, norms=%g, res_goal=%g\n", rsq,
			normsrc, rsqstop);
#endif

	if (rsq < rsqstop) { // already converged with initial vector
		*resid = rsq;
		return;
	}

	//Start iterations:
	dvecdvec = rsq;

	for (ite = 0; ite < maxit; ++ite) {
		oldrsq = rsq;
#ifdef DEBUG_TEST_PROJ
		trivialOperator(0.2,dvec,tmp);
#else
		DSq(dvec, tmp);
#endif
		//aval =|rvec_i|^2/(<dvec_i z> + sigma <dvec_i dvec_i>)
		sumtmp = 0.0;
		FORALLSITES(i, s)
		{
			ctmp = TF_dot(&(dvec[i]), &(tmp[i]));
			sumtmp += ctmp.real;
		}
		g_doublesum(&sumtmp);
		dvecz = sumtmp;
		aval = rsq / (dvecz + shift * dvecdvec);

		sumtmp = 0.0;
		FORALLSITES(i, s)
		{
			scalar_mult_sum_TF(&(dvec[i]), aval, &(dest[i])); // dest+=aval*dvec
			scalar_mult_sum_TF(&(tmp[i]), -aval, &(rvec[i])); // rvec-=aval*(tmp+shift*dvec)
			scalar_mult_sum_TF(&(dvec[i]), -aval * shift, &(rvec[i]));
			sumtmp += (double) magsq_TF(&(rvec[i]));
		}
		g_doublesum(&sumtmp);
		rsq = sumtmp;

#ifdef DEBUG_CG
		node0_printf("Mode CG it=%d, resid=%g, res_goal=%g\n", ite + 1, rsq,
				rsqstop);
#endif

		if (rsq < rsqstop) {

			// This is to ensure the correctness of the inversion.
			// If inversion is not correct this forces a complete restart of the iterations.
#ifdef DEBUG_TEST_PROJ
			trivialOperator(0.2,dest,tmp);
#else
			DSq(dest, tmp);
#endif
			sumtmp = 0.0;
			FORALLSITES(i, s)
			{
				scalar_mult_mult_add_TF(&(tmp[i]), -1.0, &(dest[i]), -shift,
						&(rvec[i]));
				scalar_mult_sum_TF(&(src[i]), 1.0, &(rvec[i]));
				copy_TF(&(rvec[i]), &(dvec[i]));
				sumtmp += (double) magsq_TF(&(rvec[i]));
			}
			g_doublesum(&sumtmp);
			rsq = sumtmp;

#ifdef DEBUG_CG
			node0_printf("Checking Final Mode CG res=%g\n", rsq);
#endif

			if (rsq < rsqstop) {
#ifdef DEBUG_CG
				node0_printf("Mode CG it=%d, resid=%g, res_goal=%g\n", ite + 1,
						rsq, rsqstop);
#endif
				*iterations = ite + 1;
				*resid = rsq;
				return;
			}

			// complete re-initialization
			dvecdvec = rsq;
		} else {
			beta = rsq / oldrsq;

			sumtmp = 0.0;
			FORALLSITES(i, s)
			{
				scalar_mult_add_TF(&(rvec[i]), &(dvec[i]), beta, &(dvec[i]));
				sumtmp += (double) magsq_TF(&(dvec[i]));
			}
			g_doublesum(&sumtmp);
			dvecdvec = sumtmp;
		}
	}

#ifdef DEBUG_CG
	node0_printf("Mode CG Failed convergence it=%d, resid=%g, res_goal=%g\n",
			ite + 1, rsq, rsqstop);
#endif

	*iterations = ite + 1;
	*resid = rsq;
	return;

}

/// z(x)=2x/(max-min)-(max+min)/(max-min); x=(1-2 sigma^2/(DDdag+sigma^2))
/// x^2=1-4 sigma^2 Inv + 4 sigma^4 Inv Inv
// Xoperator is z(x^2)
void applyXoperator(const Real epsilon, const Real sigma, const Real residgoal,
		const unsigned long maxit, Twist_Fermion *src, Twist_Fermion *dest,
		Twist_Fermion** workv) {
	Real sig2 = sigma * sigma;
	Real sig4 = sig2 * sig2;
	Real rescale = (2.0 / (1.0 - epsilon));
	Real subtract = ((1.0 + epsilon) / (1.0 - epsilon));
	Real idfact = (rescale - subtract);
	Twist_Fermion* tmpv = workv[3];
	Real resid;
	int iter;
	register site* s;
	int i;
	invertShiftedCG(sig2, residgoal, maxit, src, tmpv, &resid, &iter, workv);
	invertShiftedCG(sig2, residgoal, maxit, tmpv, dest, &resid, &iter, workv);
	FORALLSITES(i, s)
	{
		scalar_mult_mult_add_TF(&(dest[i]), (rescale * 4.0 * sig4), &(tmpv[i]),
				-(rescale * 4.0 * sig2), &(dest[i]));
		scalar_mult_sum_TF(&(src[i]), idfact, &(dest[i]));
	}
	/*Out = (rescale * 4.0 * sigma4) * Out - (rescale * 4.0 * sigma2) * tmp_
	 + idfact * In;*/
}

void applyXoperatorSimple(const Real epsilon, const Real sigma,
		const Real residgoal, const unsigned long maxit, Twist_Fermion *src,
		Twist_Fermion *dest, Twist_Fermion** workv) {
	Real sig2 = sigma * sigma;
	Twist_Fermion* tmpv = workv[3];
	Real resid;
	int iter;
	register site* s;
	int i;
	invertShiftedCG(sig2, residgoal, maxit, src, tmpv, &resid, &iter, workv);
	FORALLSITES(i, s)
	{
		scalar_mult_add_TF(&(src[i]), &(tmpv[i]), -(2.0 * sig2), &(dest[i]));
	}
	/*			Out = In - (2.0 * sigma2) * tmp_;*/
}



// P(z(x^2))
void clensh(const unsigned int order, const Real* coeff, const Real epsilon,
		const Real sigma, const Real residgoal, const unsigned long maxit,
		Twist_Fermion *src, Twist_Fermion *dest, Twist_Fermion** workv) {
	register site* s;
	int i;
	unsigned int k;
	Twist_Fermion* bp2;
	Twist_Fermion* bp1;
	Twist_Fermion* bn;
	Twist_Fermion* tmp;

	/// This uses the Clenshaw-Algorithm: \f$ p_n(D) v=sum_n a_n D^n v=(a_0+D(a_1+D(a_2+\ldots(a_{n-2}+D(a_{n-1} +a_{n}D)))))v \f$
	/// or:\f$ v_0=a_n v;\; v_1=Dv_0+a_{n-1}v;\; v_2=Dv_1+a_{n-2}v;\;\ldots\; v_n=p_n(D)v \f$

	if (order == 0) {
		FORALLSITES(i, s)
		{
			clear_TF(&(dest[i]));
		}
		return;
	}

	if (order == 1) {
		FORALLSITES(i, s)
		{
			scalar_mult_TF(&(src[i]), coeff[0], &(dest[i]));
		}
		return;
	}
#ifdef DEBUG_TEST_CLENSH
	trivialOperator(0.4,src,dest);
#else
	applyXoperator(epsilon, sigma, residgoal, maxit, src, dest, workv);
#endif

	if (order == 2) {
		FORALLSITES(i, s)
		{
			scalar_mult_mult_add_TF(&(src[i]), coeff[0], &(dest[i]), coeff[1],
					&(dest[i]));
		}
		return;
	}
	bp2 = dest;
	bp1 = workv[4];
	bn = workv[5];
	FORALLSITES(i, s)
	{
		scalar_mult_mult_add_TF(&(src[i]), coeff[order - 2], &(bp2[i]),
				2.0 * coeff[order - 1], &(bp2[i]));
	}
#ifdef DEBUG_TEST_CLENSH
	trivialOperator(0.4,bp2,bp1);
#else
	applyXoperator(epsilon, sigma, residgoal, maxit, bp2, bp1, workv);
#endif

	if (order == 3) {
		FORALLSITES(i, s)
		{
			scalar_mult_mult_add_TF(&(src[i]),
					coeff[order - 3] - coeff[order - 1], &(bp1[i]), 1.0,
					&(dest[i]));
		}
		return;
	}

	FORALLSITES(i, s)
	{
		scalar_mult_mult_add_TF(&(src[i]), coeff[order - 3] - coeff[order - 1],
				&(bp1[i]), 2.0, &(bp1[i]));
	}

	if (order == 4) {
#ifdef DEBUG_TEST_CLENSH
		trivialOperator(0.4,bp1,bn);
#else
		applyXoperator(epsilon, sigma, residgoal, maxit, bp1, bn, workv);
#endif
		FORALLSITES(i, s)
		{
			scalar_mult_sum_TF(&(src[i]), coeff[0], &(bn[i]));
			scalar_mult_add_TF(&(bn[i]), &(bp2[i]), -1.0, &(dest[i]));
		}
		return;
	}

	for (k = order - 4; k--;) {
#ifdef DEBUG_TEST_CLENSH
		trivialOperator(0.4,bp1,bn);
#else
		applyXoperator(epsilon, sigma, residgoal, maxit, bp1, bn, workv);
#endif
		FORALLSITES(i, s)
		{
			scalar_mult_mult_add_TF(&(src[i]), coeff[k + 1], &(bn[i]), 2.0,
					&(bn[i]));
			scalar_mult_sum_TF(&(bp2[i]), -1.0, &(bn[i]));
		}
		tmp = bp2;
		bp2 = bp1;
		bp1 = bn;
		bn = tmp;
	}
#ifdef DEBUG_TEST_CLENSH
	trivialOperator(0.4,bp1,bn);
#else
	applyXoperator(epsilon, sigma, residgoal, maxit, bp1, bn, workv);
#endif

	FORALLSITES(i, s)
	{
		scalar_mult_sum_TF(&(src[i]), coeff[0], &(bn[i]));
		scalar_mult_add_TF(&(bn[i]), &(bp2[i]), -1.0, &(dest[i]));
	}

}

int calculateGLMethod(const unsigned int noshifts, const Real* shiftparam,
		const Real rescalefact, Real* result, Real* staterr,
		const unsigned int nest,
		const unsigned int order, const Real* coeff, const Real epsilon,
		const Real residgoal, const unsigned int maxit) {
	register site* s;
	int i;
	int k, l;
	Real rescsig=1.0;
	Real sumtmp=0.0;
	Real sqtmp=0.0;
	Twist_Fermion *v0 = malloc(sites_on_node * sizeof(Twist_Fermion));
	Twist_Fermion *v1 = malloc(sites_on_node * sizeof(Twist_Fermion));
	Twist_Fermion *v2 = malloc(sites_on_node * sizeof(Twist_Fermion));
	Twist_Fermion * workv[6];
	workv[0] = malloc(sites_on_node * sizeof(Twist_Fermion));
	workv[1] = malloc(sites_on_node * sizeof(Twist_Fermion));
	workv[2] = malloc(sites_on_node * sizeof(Twist_Fermion));
	workv[3] = malloc(sites_on_node * sizeof(Twist_Fermion));
	workv[4] = malloc(sites_on_node * sizeof(Twist_Fermion));
	workv[5] = malloc(sites_on_node * sizeof(Twist_Fermion));
#ifdef DEBUG_TEST_CLENSH
	rand_TFsourceZ2(v0);
	clensh(order, coeff, epsilon, rescsig, residgoal, maxit, v0, v1,
						workv);
	sumtmp = 0.0;
	FORALLSITES(i, s)
	{
		sumtmp += (double) magsq_TF(&(v1[i]));
	}
	sumtmp*=(1.0 / ((double) DIMF * (1 + 5 + NPLAQ) * sites_on_node));
	g_doublesum(&sumtmp);
	node0_printf("TESTMODE Clensh: P^2(z)= %lg\n",sumtmp);
	rescsig=0.2;
	applyXoperatorSimple(epsilon,rescsig,residgoal,maxit,v0,v1,workv);
	sumtmp = 0.0;
	FORALLSITES(i, s)
	{
		sumtmp += (double) magsq_TF(&(v1[i]));
	}
	sumtmp*=(1.0 / ((double) DIMF * (1 + 5 + NPLAQ) * sites_on_node));
	g_doublesum(&sumtmp);
	node0_printf("TESTMODE Clensh: xOperator^2(0.2,0.4)= %lg\n",sumtmp);
	applyXoperator(epsilon,rescsig,residgoal,maxit, v0,v1,workv);
	sumtmp = 0.0;
	FORALLSITES(i, s)
	{
		sumtmp += (double) magsq_TF(&(v1[i]));
	}
	sumtmp*=(1.0 / ((double) DIMF * (1 + 5 + NPLAQ) * sites_on_node));
	g_doublesum(&sumtmp);
	node0_printf("TESTMODE Clensh: F(z(xOperator(0.2,0.4)))^2= %lg\n",sumtmp);
	free(v0);
	free(v1);
	free(v2);
	free(workv[0]);
	free(workv[1]);
	free(workv[2]);
	free(workv[3]);
	free(workv[4]);
	free(workv[5]);
	return (0);

#endif

	for (k = 0; k < noshifts; k++) {
		rescsig = shiftparam[k] / rescalefact;
		result[k] = 0.0;
		sqtmp = 0.0;
		for (l = 0; l < nest; l++) {
			rand_TFsourceZ2(v0);

			/// h(x)=0.5(1-xP(z(x^2)));
			/// z(x)=2x/(max-min)-(max+min)/(max-min);
			/// x=(1-2sigma^2/(DDdag+sigma^2))
			/// h(x) applied
			clensh(order, coeff, epsilon, rescsig, residgoal, maxit, v0, v1,
					workv);
			applyXoperatorSimple(epsilon, rescsig, residgoal, maxit, v1, v2,
					workv);
			FORALLSITES(i, s)
			{
				scalar_mult_mult_add_TF(&(v0[i]), 0.5, &(v2[i]), -0.5,
						&(v1[i]));
			}
			clensh(order, coeff, epsilon, rescsig, residgoal, maxit, v1, v0,
					workv);
			applyXoperatorSimple(epsilon, rescsig, residgoal, maxit, v0, v2,
					workv);
			sumtmp = 0.0;
			FORALLSITES(i, s)
			{
				scalar_mult_mult_add_TF(&(v1[i]), 0.5, &(v2[i]), -0.5,
						&(v0[i]));
				sumtmp += (double) magsq_TF(&(v0[i]));
			}
			sumtmp*=(1.0 / ((double) DIMF * (1 + 5 + NPLAQ) * sites_on_node));
			g_doublesum(&sumtmp);
			result[k] += sumtmp / ((Real) nest);
			sqtmp += sumtmp * sumtmp / ((Real) nest);

		}
		staterr[k] = sqrt((sqtmp - result[k] * result[k])*((Real)nest)/((Real) nest-1.0));
		node0_printf("Finished GL ModeN %d %g %g %g %g\n", k, shiftparam[k],
				result[k], staterr[k], rescsig);
	}
	free(v0);
	free(v1);
	free(v2);
	free(workv[0]);
	free(workv[1]);
	free(workv[2]);
	free(workv[3]);
	free(workv[4]);
	free(workv[5]);
	return (0);
}

// -----------------------------------------------------------------
//
int calculate_coeff(unsigned int Nest, unsigned int order, Real lambda_min,
		Real lambda_max, Real** ckcoeff) {
	register site* s;
	int i, j, k;
	Twist_Fermion *randsrc = malloc(sites_on_node * sizeof(Twist_Fermion));
	Twist_Fermion * v0 = malloc(sites_on_node * sizeof(Twist_Fermion));
	Twist_Fermion * vk = malloc(sites_on_node * sizeof(Twist_Fermion));
	Twist_Fermion * vkp1 = malloc(sites_on_node * sizeof(Twist_Fermion));
	Twist_Fermion * vkm1 = malloc(sites_on_node * sizeof(Twist_Fermion));
	Twist_Fermion * tmp = malloc(sites_on_node * sizeof(Twist_Fermion));
	Real diff = lambda_max - lambda_min;
	Real diffp = lambda_max + lambda_min;
	Real fakt1 = 2.0 / diff;
	Real fakt2 = -diffp / diff;
	Real vnorm =
			(1.0 / ((double) DIMF * (1 + 5 + NPLAQ)  * sites_on_node));
	Real vnorm2=1.0/((double)Nest);
	Real sumtmp = 0.0;
	complex ctmp;

	if (*ckcoeff == NULL) {
		*ckcoeff = malloc(order * sizeof(Real));
	}

	// Check memory allocations
	if (!randsrc || !v0 || !vk || !vkp1 || !vkm1 || !tmp || !*ckcoeff) {
		node0_printf("ERROR: Allocation in calculate_coef failed!\n");
		exit(1);
	}

	for (j = 0; j < order; j++) {
		(*ckcoeff)[j] = 0.0;
	}

	for (k = 0; k < Nest; k++) {
		rand_TFsourceZ2(randsrc);
		sumtmp = 0.0;
		FORALLSITES(i, s)
		{
			copy_TF(&(randsrc[i]), &(v0[i]));
			sumtmp += (double) magsq_TF(&(v0[i]));
		}
		sumtmp = vnorm * sumtmp;
		g_doublesum(&sumtmp);
		(*ckcoeff)[0] += sumtmp*vnorm2;

		if (order > 1) {
			DSq(v0, tmp);
			sumtmp = 0.0;
			FORALLSITES(i, s)
			{
				scalar_mult_mult_add_TF(&(tmp[i]), fakt1, &(v0[i]), fakt2,
						&(vk[i]));
				ctmp = TF_dot(&(v0[i]), &(vk[i]));
				sumtmp += ctmp.real;
				copy_TF(&(v0[i]), &(vkm1[i]));
			}
			sumtmp*=vnorm;
			g_doublesum(&sumtmp);
			(*ckcoeff)[1] += vnorm2 * sumtmp;
			for (j = 2; j < order; j++) {
				DSq(vk, tmp);
				sumtmp = 0.0;
				FORALLSITES(i, s)
				{
					scalar_mult_mult_add_TF(&(tmp[i]), 2.0 * fakt1, &(vk[i]),
							2.0 * fakt2, &(vkp1[i]));
					scalar_mult_sum_TF(&(vkm1[i]), -1.0, &(vkp1[i]));
					ctmp = TF_dot(&(v0[i]), &(vkp1[i]));
					sumtmp += ctmp.real;
					copy_TF(&(vk[i]), &(vkm1[i]));
					copy_TF(&(vkp1[i]), &(vk[i]));
				}
				sumtmp*=vnorm;
				g_doublesum(&sumtmp);
				(*ckcoeff)[j] += vnorm2 * sumtmp;
			}
		} /* if order >1 */
		node0_printf("Finished Estimator %d\n", k);
	} /* end loop over estimators */
	/*for (j = 0; j < order; j++) {
		node0_printf("ChebCoeff: %.8g\n", (*ckcoeff)[j]);
	}*/
	/*free(ckcoeff); has to be called outside*/
	free(randsrc);
	free(v0);
	free(vk);
	free(vkp1);
	free(vkm1);
	free(tmp);
	return (1);
}

// -----------------------------------------------------------------
