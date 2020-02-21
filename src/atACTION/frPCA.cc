#include <ATACtion.h>

namespace ATACtion {
	mat randNorm(int l, int m, int seed) {
		std::default_random_engine gen (seed);	
		std::normal_distribution<double> unif(0.0, 1.0);
			
		mat R(l, m);
		for (register int j = 0; j < m; j++) {
			for(register int i = 0; i < l; i++) {
				R(i, j) = unif(gen);
			}
		}
		return R;
	}
	
	
	field<mat> eigSVD(mat A) {
		int n = A.n_cols;
		mat B = trans(A)*A;
				
		vec d;
		mat V;
		eig_sym( d, V, B );		
		d = sqrt(d);
		
		// Compute U
		sp_mat S(n, n);
		S.diag() = 1 / d; 
		mat U = (S*trans(V))*trans(A);
		U = trans(U);
		
		field<mat> out(3);
		
		out(0) = U;
		out(1) = d;
		out(2) = V;
		
		return(out);
	}


/*
	field<mat> eigSVD(sp_mat A, int k) {
		int n = A.n_cols
		sp_mat B = trans(A)*A;
		
		vec d;
		mat V;
		eigs_sym( d, V, B, k) 
		d = sqrt(d)
		
		// Compute U
		sp_mat S(k, k);
		S.diag() = 1 / d; 
		mat U = (S*trans(V))*trans(A);
		U = trans(U);
		
		field out(3);
		
		out(0) = U;
		out(1) = d;
		out(2) = V;
	}
*/

	// From: Xu Feng, Yuyang Xie, and Yaohang Li, "Fast Randomzied PCA for Sparse Data," in Proc. the 10th Asian Conference on Machine Learning (ACML), Beijing, China, Nov. 2018.
	field<mat> frSVD(sp_mat &A, int dim, int iters, int seed = 0) {	
		int s = 5;

		int m = A.n_rows;
		int n = A.n_cols;
		
		printf("\t\tRunning randomized PCA. Matrix size: %d x %d (# iters = %d)\n", m, n, iters); fflush(stdout);
				
		vec S;
		mat Q, L, U, V;
		field<mat> SVD_out;
		
		if(m < n) {
			printf("\t\t\tInitializing PCA (mode 1) ... ");
			Q = randNorm(n, dim+s, seed);
			Q = A*Q;
			if (iters == 0) {
				SVD_out = eigSVD(Q);
				Q = SVD_out(0);					
			}
			else {
				lu(L, U, Q);
				Q = L;
			}
			printf("done\n");
			
			for (int i = 1; i <= iters; i++) {
				printf("\t\t\t\tIter %d/%d ... ", i, iters);				
				
				if (i == iters) {
					SVD_out = eigSVD(A*(trans(A)*Q));
					Q = SVD_out(0);									
				}
				else {
					lu(L, U, A*(trans(A)*Q));
					Q = L;
				}
				printf("done\n");				
			}
			
			SVD_out = eigSVD(trans(A)*Q);
			V = SVD_out(0);
			S = vec(SVD_out(1));
			U = SVD_out(2);
			
			U = Q*fliplr(U.cols(s, dim+s-1));
			V = fliplr(V.cols(s, dim+s-1));
			S = flipud(S(span(s, dim+s-1)));
		}
		else {
			printf("\t\t\tInitializing PCA (mode 2) ... ");				
			Q = randNorm(m, dim+s, seed);
			Q = trans(A)*Q;
			if (iters == 0) {
				SVD_out = eigSVD(Q);
				Q = SVD_out(0);					
			}
			else {
				lu(L, U, Q);
				Q = L;
			}
			printf("done\n");
			
			for (int i = 1; i <= iters; i++) {
				printf("\t\t\t\tIter %d/%d ... ", i, iters);				
				
				if (i == iters) {
					SVD_out = eigSVD(trans(A)*(A*Q));
					Q = SVD_out(0);									
				}
				else {
					lu(L, U, trans(A)*(A*Q));
					Q = L;
				}
				printf("done\n");				
				
			}
			
			SVD_out = eigSVD(A*Q);
			U = SVD_out(0);
			S = vec(SVD_out(1));
			V = SVD_out(2);
						
			
			U = fliplr(U.cols(s, dim+s-1));
			V = Q*fliplr(V.cols(s, dim+s-1));
			S = flipud(S(span(s, dim+s-1)));
		}		
	
	
		field<mat> out(3);		
		out(0) = U;
		out(1) = S;
		out(2) = V;
		
		printf("\t\t\tdone\n");	fflush(stdout);			
		
		return(out);
	}


	void gram_schmidt(mat& A) {
		for(uword i = 0; i < A.n_cols; ++i) {
			for(uword j = 0; j < i; ++j) {
				double r = dot(A.col(i), A.col(j));
				A.col(i) -= r * A.col(j);
			}
			
			double col_norm = norm(A.col(i), 2);
			
			if(col_norm < 1E-4) {
				for(uword k = i; k < A.n_cols; ++k)
					A.col(k).zeros();

				return;
			}
			A.col(i) /= col_norm;
		}
	}

	Projection reducedKernel(sp_mat &profile, int PCA_dim, int iter = 3, int seed = 0) {			
		int n = profile.n_rows;
		profile = normalise(profile, 2);    

		Projection projection;		

		printf("\tRunning ACTION+PCA. Matrix size: %d x %d\n", profile.n_rows, profile.n_cols); fflush(stdout);

		mat A, B;
		
		// Update 1: Orthogonalize columns w.r.t. background (mean)
		vec mu = vec(mean(profile, 1));
		vec v = mu / norm(mu, 2);
		vec a1 = v;
		vec b1 = -trans(profile)*v;
		
		// Update 2: Center columns of orthogonalized matrix before performing SVD (i.e., PCA)
		vec c = vec(trans(mean(profile, 0)));
		double a1_mean = mean(a1);
		vec a2 = ones(profile.n_rows);
		vec b2 = -(a1_mean*b1 + c);

		A = join_rows(a1, a2);
		B = join_rows(b1, b2);

		
		printf("\tPerform SVD on the original matrix\n"); fflush(stdout);
		field<mat> SVD_results = frSVD(profile, PCA_dim, iter, seed);	

		mat U = SVD_results(0);
		vec s = SVD_results(1);
		mat V = SVD_results(2);
		

//		printf("U: %d x %d, V = %d x %d, s = %d\n", U.n_rows, U.n_cols, V.n_rows, V.n_cols, s.n_elem);
		printf("\tUpdate SVD ..."); fflush(stdout);
		vec s_prime;
		mat U_prime, V_prime;
			
		mat M = U.t()*A; 
		mat A_ortho_proj = A - U*M;   
		mat P = A_ortho_proj;// = orth(A_ortho_proj);
		gram_schmidt(P);		
		mat R_P = P.t()*A_ortho_proj;
		
		
		mat N = V.t()*B; 
		mat B_ortho_proj = B - V*N; 
		mat Q = B_ortho_proj; //orth(B_ortho_proj); 
		gram_schmidt(Q);
		mat R_Q = Q.t()*B_ortho_proj;	
		
		mat K1 = zeros(s.n_elem+A.n_cols, s.n_elem+A.n_cols);
		for(int i = 0; i < s.n_elem; i++) {
			K1(i, i) = s(i);
		}

		mat K2 = join_vert(M, R_P)*trans(join_vert(N, R_Q));

		
		mat K = K1 + K2;

		svd( U_prime, s_prime, V_prime, K );
		
		mat U_updated = join_horiz(U, P)*U_prime;
		mat V_updated = join_horiz(V, Q)*V_prime;
		printf("done.\n"); fflush(stdout);
		
		
		projection.S_r = trans(V_updated.cols(0, PCA_dim-1));
		for(int i = 0; i < PCA_dim; i++) {
			projection.S_r.row(i) *= s_prime(i);
		}
		projection.V = U_updated.cols(0, PCA_dim-1);	
		
		vec sigma_sq = square(s_prime);

		projection.lambda = sigma_sq / (n-1);

		vec explained_var = cumsum(sigma_sq) / sum(sigma_sq);
		projection.exp_var = explained_var(span(0, PCA_dim-1));
		
		return projection;	
	}
}
