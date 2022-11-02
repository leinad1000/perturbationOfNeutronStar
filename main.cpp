#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

void initialValues(double& p,
                   double& m, 
	               double& r, 
	               double& phi, 
	               double& rho, 
	               double& dphi, 
	               double pStart, 
	               double rStart, 
	               double phiStart, 
	               double rho_c, 
	               double dphiStart);

void finding_n(int& n,
	           double& p,
	           double& m,
	           double& r,
	           double& phi, 
	           double& dphi, 
	           double& rho, 
	           double pStart, 
	           double toleranse,
	           double mu, 
	           double lambda, 
	           double M,
	           double rStep,
	           double phiVacuum,
	           int stepWhenPhiStrikes,
	           double faktor, 
	           double xKonstant,
	           double pKonstant,
	           double epsilon);

void equilibriumValues(int n, 
                       double& p, 
                       double& m, 
                       double& r, 
                       double& phi, 
                       double& dphi, 
                       double& rho, 
                       double pStart, 
                       double toleranse, 
                       double mu, 
                       double lambda, 
                       double M, 
                       double rStep,
                       double phiVacuum, 
                       int stepWhenPhiStrikes, 
                       double faktor, 
                       double* rVector, 
                       double* pVector, 
                       double* dpVector, 
                       double* ddpVector, 
                       double* mVector, 
                       double* dmVector, 
                       double* rhoVector, 
                       double* drhoVector, 
                       double* phiVector,
                       double* dphiVector,
                       double* ddphiVector, 
                       double* e2betaVector, 
                       double* dbetaVector, 
                       double* dalphaVector, 
                       double* ddalphaVector,
                       double* gVector, 
                       double* fVector, 
                       double* dVvector, 
                       double* ddVvector, 
                       double* Vvector, 
                       double* kVector, 
                       double* dkVector,
                       double* ddrhoVector,
                       double epsilon, 
                       double A, 
                       double B,
                       double xKonstant,
                       double pKonstant);

void findingAlpha(int n, 
                  double* dphiVector,
                  double* ddphiVector,
                  double* pVector, 
                  double* dpVector,
                  double* rhoVector,
                  double* dVvector, 
                  double* rVector, 
                  double* Vvector,
                  double* mVector,
                  double* alphaVector);

void functions(int n, 
               double* rVector,
               double* pVector,
               double* dpVector, 
               double* ddpVector,
               double* mVector, 
               double* dmVector, 
               double* rhoVector,
               double* drhoVector,
               double* phiVector, 
               double* dphiVector,
               double* ddphiVector, 
               double* e2betaVector,
               double* dbetaVector,
               double* alphaVector,
               double* dalphaVector,
               double* ddalphaVector,  
               double* gVector, 
               double* fVector,
               double* dVvector, 
               double* ddVvector,
               double* Vvector,
               double* a,
               double* c,
               double* d, 
               double* aTilde,
               double* cTilde, 
               double* dTilde, 
               double* omegaTerm,
               double* bTerm, 
               double* bBarTerm, 
               double* nevner,
               double* ddrhoVector,
               double* kVector,
               double* dkVector,
               double B);

void eigenvalues(double stepWhenPhiStrikes,
                 double rStep,
                 double xiMax,
                 int n,
                 int N,
                 double* a,
                 double* b,
                 double* c,
                 double* d,
                 double* aTilde,
                 double* bTilde,
                 double* cTilde,
                 double* dTilde,
                 double* omegaTerm,
                 double* bTerm,
                 double* bBarTerm,
                 double* xi,
                 double* deltaPhi,
                 double* xiZero,
                 double* omega2,
                 double* nevner);

void tellerEgenverdier(int N, 
                       int& antallEgenverdier, 
                       double* xiZero);

void finnerEgenverdier(double stepWhenPhiStrikes,
                       double xiMax,
                       double rStep,
                       int antallEgenverdier,
                       int N,
                       int n,
                       double* xiZero,
                       double* omega2,
                       double* a,
                       double* b,
                       double* c,
                       double* d,
                       double* bTerm,
                       double* aTilde,
                       double* bTilde,
                       double* cTilde,
                       double* dTilde,
                       double* bBarTerm,
                       double* omegaTerm,
                       double* nevner,
                       double* xi,
                       double* deltaPhi,
                       double* egenverdiVektor,
                       double* besteXiNullVektor);

int main() {

	// Åpner ut-fil
	ofstream myfile;
	myfile.open("TestFil.txt");

	// Variabler   
	double 
		sum1,
		sum2,
		pi,
		toleranse,
		rho_c,
		rStart,
		rStep,
		xiMax,
		pStart,
		phiStart,
		dphiStart,
		phiVacuum,
		mu,
		lambda,
		p,
		m,
		r,
		phi,
		rho,
		dphi,
		M,
		H0,
		mNeutron,
		hBar,
		step_faktor,
		epsilon,
		faktorMin,
		faktorMax,
		omegaMin,
		omegaMax,
		omegaStep,
		faktor,
		x,
		xKonstant,
		pKonstant,
		A,
		B,
		rho_cMin,
		rho_cMax,
		rho_cStep;

	int N, 
		n, 
		n1, 
		n2, 
		stepWhenPhiStrikes, 
		antallEgenverdier, 
		nn;

	// Konstanter til simuleringen
	rStart = 6.770190656e-7; // 1 mm
    rStep = 6.770190656e-4; // 1 m
	xiMax = 6.770190656e-4; // Perturbasjonen
	toleranse = 1e-13; // Små tettheter (ved overflaten til stjerna)
	epsilon = 1e-13; // Parameter for å avslutte løkke
	nn = 1e2; // Antall datapunkter

	// Naturlige konstanter, i enheter hvor c=G=M_sol=1
	H0 = (74.0/69.8) * 1.11450822e-23; // Hubbles konstant, når den er målt til å være 69.8 km/s/Mpc
	mNeutron = 0.842095273e-57; // nøytronmassen
	hBar = 1.19735126e-76; // Plancks reduserte konstant
	pi = acos(-1.0L);

	// Samle-konstanter
	pKonstant = pow(mNeutron, 4.0) / (8 * pow(pi, 2.0) * pow(hBar, 3.0));
	xKonstant = pow(3 * pow(pi,2.0) * pow(hBar,3.0) / (pow(mNeutron, 4.0)), (1.0 / 3.0));
	A = pow(mNeutron, (8.0 / 3.0)) / (8.0 * pow(3.0, (2.0 / 3.0)) * pow(hBar, 2.0) * pow(pi, (4.0 / 3.0)));
	B = pow(pow(pi, 2.0) * pow(hBar, 3.0) / (9.0 * pow(mNeutron, 4.0)), (1.0 / 3.0));

	// Startverdier
	dphiStart = 0.0;

	// Sentral-tetthet
	// rho_cMin = 9.7198e-5; // Minste verdi
	// rho_cMax = 4.8599e-2; // Største verdi
    rho_cMin = 0.008425790262626263; // Minste verdi
    rho_cMax = 0.009405624646464647; // Største verdi
	rho_cStep = (rho_cMax - rho_cMin) / ((double)(nn - 1));

	// Symmetron-konstanter
	M = (1e-4) / (sqrt(8 * pi));
	mu = H0 / (sqrt(16.0 * pi) * M);
	lambda = pow(H0, 2.0) / (64.0 * pow(pi, 2.0) * pow(M, 6));
	phiVacuum = mu / sqrt(lambda);
    
    
    
	// Simulering
	for (int ii = 0; ii < nn; ii++) {
        
		// Nullstiller
		stepWhenPhiStrikes = 1e6;
		n = 0;
		N = 3;
		faktor = 1.0;

		// Ny sentral-tetthet for hvert steg i simuleringen, med tilhørende sentral-trykk
		rho_c = rho_cMin + ii * rho_cStep;
        x = xKonstant * pow(rho_c, 1.0 / 3.0);
		pStart = pKonstant*(x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));

        // Finner radien på stjerna uten skalarfelt tilstede
		phiStart = 0;
		initialValues(p, m, r, phi, rho, dphi, pStart, rStart, phiStart, rho_c, dphiStart);
		finding_n(n, p, m, r, phi, dphi, rho, pStart, toleranse, mu, lambda, M, rStep, phiVacuum, stepWhenPhiStrikes, faktor, xKonstant, pKonstant, epsilon);
        n1 = n;
        
        // Finner radien på stjerna med skalarfelt tilstede allerede fra første steg
        phiStart = phiVacuum * (1e-250);
        n = 0;
        initialValues(p, m, r, phi, rho, dphi, pStart, rStart, phiStart, rho_c, dphiStart);
        finding_n(n, p, m, r, phi, dphi, rho, pStart, toleranse, mu, lambda, M, rStep, phiVacuum, stepWhenPhiStrikes, faktor, xKonstant, pKonstant, epsilon);
        n2 = n;

        // Steget der hvor symmetronet skal settes forskjellig fra null
        stepWhenPhiStrikes = n1 - n2;
       
        phiStart = 0.0; // Setter phi i sentrum av stjerna til å være lik null igjen i det følgende

        // Finner ut hva phi skal settes til når den settes forskjellig fra null (ved stepWhenPhiStrikes) for at phi ved overflaten skal være lik phiVacuum.
        double* faktorVector = new double[N];
        double* phiEndVector = new double[N];

        faktorMin = 0.0;
        faktorMax = 1.0;

        while (fabs(1.0 - faktorMin / faktorMax) > epsilon) {
            step_faktor = (faktorMax - faktorMin) / ((double)(N - 1));
            for (int i = 0; i < N; i++) {
                faktorVector[i] = faktorMin + i * step_faktor;
                faktor = faktorVector[i];
                n = 0; // Nullstiller
                initialValues(p, m, r, phi, rho, dphi, pStart, rStart, phiStart, rho_c, dphiStart); // Nullstiller
                finding_n(n, p, m, r, phi, dphi, rho, pStart, toleranse, mu, lambda, M, rStep, phiVacuum, stepWhenPhiStrikes, faktor, xKonstant, pKonstant, epsilon);
                phiEndVector[i] = phi;
            }
            for (int i = 0; i < N; i++) {
                if (phiEndVector[i] < phiVacuum) {
                    faktorMin = faktorVector[i];
                    faktorMax = faktorVector[i + 1];
                }
            }
        } 

        faktor = 0.5 * (faktorMin + faktorMax); // faktor som forteller hva første phi forskjellig fra null skal settes til

        // Finner antall datapunkter, det vil si "n"
        n = 0; // Nullstiller
        initialValues(p, m, r, phi, rho, dphi, pStart, rStart, phiStart, rho_c, dphiStart);
        // Regner ut "n" med den korrekte startverdien av Symmetronet
        finding_n(n, p, m, r, phi, dphi, rho, pStart, toleranse, mu, lambda, M, rStep, phiVacuum, stepWhenPhiStrikes, faktor, xKonstant, pKonstant, epsilon);
    
        // Vektorer som skal inneholde likevektsverdiene
        double* rVector = new double[n];
        double* pVector = new double[n];
        double* dpVector = new double[n];
        double* ddpVector = new double[n];
        double* mVector = new double[n];
        double* dmVector = new double[n];
        double* rhoVector = new double[n];
        double* drhoVector = new double[n];
        double* ddrhoVector = new double[n];
        double* phiVector = new double[n];
        double* dphiVector = new double[n];
        double* ddphiVector = new double[n];
        double* e2betaVector = new double[n];
        double* dbetaVector = new double[n];
        double* alphaVector = new double[n];
        double* dalphaVector = new double[n];
        double* ddalphaVector = new double[n];
        double* gVector = new double[n];
        double* fVector = new double[n];
        double* dVvector = new double[n];
        double* ddVvector = new double[n];
        double* Vvector = new double[n];
        double* kVector = new double[n];
        double* dkVector = new double[n];
  
        // Nullstiller
        initialValues(p, m, r, phi, rho, dphi, pStart, rStart, phiStart, rho_c, dphiStart); 
        // Finner likevektsverdiene
        equilibriumValues(n, p, m, r, phi, dphi, rho, pStart, toleranse, mu, lambda, M, rStep, phiVacuum, stepWhenPhiStrikes, faktor, rVector, pVector, dpVector, ddpVector, mVector, dmVector, rhoVector, drhoVector, phiVector, dphiVector, ddphiVector, e2betaVector, dbetaVector, dalphaVector, ddalphaVector, gVector, fVector, dVvector, ddVvector, Vvector, kVector, dkVector, ddrhoVector, epsilon, A, B, xKonstant, pKonstant); 
        // Regner ut alpha til slutt
        findingAlpha(n, dphiVector, ddphiVector, pVector, dpVector, rhoVector, dVvector, rVector, Vvector, mVector, alphaVector); 
        
        // Nye vektorer som skal brukes til å løse de to koblede differensialligningene. Hver løsning gir en egenverdi.
        double* a = new double[n];
        double* b = new double[n];
        double* c = new double[n];
        double* d = new double[n];
        double* aTilde = new double[n];
        double* bTilde = new double[n];
        double* cTilde = new double[n];
        double* dTilde = new double[n];
        double* omegaTerm = new double[n];
        double* bTerm = new double[n];
        double* bBarTerm = new double[n];
        double* nevner = new double[n];

        // Regner ut verdiene til funksjonene definert like over
        functions(n, rVector, pVector, dpVector, ddpVector, mVector, dmVector, rhoVector, drhoVector, phiVector, dphiVector, ddphiVector, e2betaVector, dbetaVector, alphaVector, dalphaVector, ddalphaVector, gVector, fVector, dVvector, ddVvector, Vvector, a, c, d, aTilde, cTilde, dTilde, omegaTerm, bTerm, bBarTerm, nevner, ddrhoVector, kVector, dkVector, B);

        N = 1e4; // Øker N for å sørge for at jeg finner alle egenverdiene i intervallet

        // Vektorer som skal inneholde forslag til egenverdi og resulterende xi-null. Når xi-null = 0, er en egenverdi funnet.
        double* omega2 = new double[N];
        double* xiZero = new double[N];
        // Vektorer som skal inneholde perturbasjonene.
        double* xi = new double[n];
        double* deltaPhi = new double[n];

        // Grid for egenverdier
        omegaMin = -0.1; // Nedre grense
        omegaMax = 0.1; // Øvre grense
        omegaStep = (omegaMax - omegaMin) / ((double)(N - 1));
        for (int i = 0; i < N; i++) {
            omega2[i] = omegaMin + i * omegaStep;
        }

        // Sjekker om det er egenverdier til stede
        eigenvalues(stepWhenPhiStrikes, rStep, xiMax, n, N, a, b, c, d, aTilde, bTilde, cTilde, dTilde, omegaTerm, bTerm, bBarTerm, xi, deltaPhi, xiZero, omega2, nevner);

        // Teller hvor mange egenverdier som er funnet i det valgte intervallet
        tellerEgenverdier(N, antallEgenverdier, xiZero);

        // Definerer vektor som skal inneholde alle egenverdiene som er funnet i intervallet
        double* egenverdiVektor = new double[antallEgenverdier];
        double* besteXiNullVektor = new double[antallEgenverdier];

        // Finner egenverdiene innenfor ønsket presisjon, og tar vare på de
        finnerEgenverdier(stepWhenPhiStrikes, xiMax, rStep, antallEgenverdier, N, n, xiZero, omega2, a, b, c, d, bTerm, aTilde, bTilde, cTilde, dTilde, bBarTerm, omegaTerm, nevner, xi, deltaPhi, egenverdiVektor, besteXiNullVektor);

        /*
        // Skriver de fundamentale egenverdiene til fil
        myfile << setw(20) << setprecision(16) << rho_c << "\t";
        myfile << setw(20) << setprecision(16) << egenverdiVektor[0] << "\t";
        myfile << setw(20) << setprecision(16) << rVector[n - 1] << "\t";
        myfile << setw(20) << setprecision(16) << mVector[n - 1] << endl;
        */

        /*
        // Skriver Symmetron-profile til fil for en utvalgt sentral-tetthet
        for (int i = 0; i < n; i++) {
            myfile << setw(20) << setprecision(16) << rVector[i] << "\t"; // Radius i enheter hvor c=G=M_sol=1
            myfile << setw(20) << setprecision(16) << (rVector[i] / 0.6770190656) << "\t"; // Radius i antall kilometer
            myfile << setw(20) << setprecision(16) << phiVector[i] << endl; // Symmetronet inni stjerna
        }
        */

        /*
        // Skriver Masse-profilen til fil
        myfile << setw(20) << setprecision(16) << r << "\t"; // Radius i enheter hvor c=G=M_sol=1
        myfile << setw(20) << setprecision(16) << ((double)(n)/1000.0) << "\t"; // Radius i antall kilometer
        myfile << setw(20) << setprecision(16) << m << "\t"; // Masse i antall solmasser
        myfile << setw(20) << setprecision(16) << rho_c << endl; // sentral-tetthet
        */     
        
        // Sletter vektorer
        delete[] egenverdiVektor;
        delete[] besteXiNullVektor;
        delete[] omega2;
        delete[] xiZero;
        delete[] xi;
        delete[] deltaPhi;
        delete[] a;
        delete[] b;
        delete[] c;
        delete[] d;
        delete[] aTilde;
        delete[] bTilde;
        delete[] cTilde;
        delete[] dTilde;
        delete[] omegaTerm;
        delete[] bTerm;
        delete[] bBarTerm;
        delete[] nevner;
        delete[] rVector;
        delete[] pVector;
        delete[] dpVector;
        delete[] ddpVector;
        delete[] mVector;
        delete[] dmVector;
        delete[] rhoVector;
        delete[] drhoVector;
        delete[] ddrhoVector;
        delete[] phiVector;
        delete[] dphiVector;
        delete[] ddphiVector;
        delete[] e2betaVector;
        delete[] dbetaVector;
        delete[] alphaVector;
        delete[] dalphaVector;
        delete[] ddalphaVector;
        delete[] gVector;
        delete[] fVector;
        delete[] dVvector;
        delete[] ddVvector;
        delete[] Vvector;
        delete[] kVector;
        delete[] dkVector;
        delete[] faktorVector;
        delete[] phiEndVector;
    
        cout << "Antall linjer skrevet til fil er: " << ii + 1 << endl;
	}
    	
	// Stenger ut-fil
	myfile.close();
	return 0;
}

void initialValues(double& p, 
	               double& m,
	               double& r, 
	               double& phi, 
	               double& rho, 
	               double& dphi, 
	               double pStart, 
	               double rStart, 
	               double phiStart, 
	               double rho_c, 
	               double dphiStart) {
	p = pStart;
	m = 0.0;
	r = rStart;
	phi = phiStart;
	rho = rho_c;
	dphi = dphiStart;
}

void finding_n(int& n,
               double& p,
               double& m,
               double& r,
               double& phi,
               double& dphi,
               double& rho,
               double pStart,
               double toleranse,
               double mu, 
               double lambda,
               double M,
               double rStep,
               double phiVacuum, 
               int stepWhenPhiStrikes, 
               double faktor, 
               double xKonstant,
               double pKonstant,
               double epsilon) {

    double 
        y,
        dy,
        dp,
        dm,
        I,
        dI,
        mu2,
        mu4,
        M2, 
        phiAnnen, 
        phiTredje,
        phiFjerde,
        r2, 
        r3,
        pi,
        constant1,
        previousPreviousI,
        previousI, 
        previous_p, 
        previous_m,
        previous_phi,
        previous_rho, 
        previous_r,
        previous_y, 
        f1, 
        f2,
        y1,
        y2,
        y3,
        y4,
        p1,
        p2,
        p3,
        p4,
        m1,
        m2,
        m3,
        m4,
        phi1,
        phi2,
        phi3,
        phi4,
        rStepHalve, 
        rStep6,
        yAnnen,
        x, 
        rhoMax,
        rhoMin,
        rhoHalve,
        pMax, 
        pMin, 
        pHalve;

    // Konstanter som endres underveis
    phiAnnen = phi * phi;
    phiTredje = phi * phiAnnen;
    phiFjerde = phiAnnen * phiAnnen;
    r2 = r * r;
    r3 = r * r2;
    y = dphi;
    yAnnen = y * y;
    // Konstanter som ikke endres underveis
    rStepHalve = rStep / 2.0;
    rStep6 = rStep / 6.0;
    mu2 = mu * mu;
    mu4 = mu2 * mu2;
    M2 = M * M;  
    pi = acos(-1.0L);
    constant1 = mu4 / (4.0 * sqrt(lambda));
   
 

    for (int i = 0; i < 2; i++) {
        if (i == 0) {
            previous_r = r;
            previous_rho = rho;
            previous_phi = phi;
            previous_p = p;
            previous_y = y;
            previous_m = m;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            f1 = dI;
            previousPreviousI = dI * r;
            I = previousPreviousI;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p1 = dp;
            m1 = dm;
            y1 = dy;
            phi1 = y;
            p = previous_p + rStepHalve * p1;
            m = previous_m + rStepHalve * m1;
            y = previous_y + rStepHalve * y1;
            phi = previous_phi + rStepHalve * phi1;

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;
            
            r = previous_r + rStepHalve;
            r2 = r * r;
            r3 = r * r2;
            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            I = previousPreviousI + 0.5 * (dI + f1) * rStepHalve;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p2 = dp;
            m2 = dm;
            y2 = dy;
            phi2 = y;
            p = previous_p + rStepHalve * p2;
            m = previous_m + rStepHalve * m2;
            y = previous_y + rStepHalve * y2;
            phi = previous_phi + rStepHalve * phi2;

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            I = previousPreviousI + 0.5 * (dI + f1) * rStepHalve;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p3 = dp;
            m3 = dm;
            y3 = dy;
            phi3 = y;
            p = previous_p + rStep * p3;
            m = previous_m + rStep * m3;
            y = previous_y + rStep * y3;
            phi = previous_phi + rStep * phi3;

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

            r = previous_r + rStep;
            r2 = r * r;
            r3 = r * r2;
            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            I = previousPreviousI + 0.5 * (dI + f1) * rStep;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p4 = dp;
            m4 = dm;
            y4 = dy;
            phi4 = y;
            p = previous_p + rStep6 * (p1 + 2.0 * p2 + 2.0 * p3 + p4);
            m = previous_m + rStep6 * (m1 + 2.0 * m2 + 2.0 * m3 + m4);
            y = previous_y + rStep6 * (y1 + 2.0 * y2 + 2.0 * y3 + y4);
            phi = previous_phi + rStep6 * (phi1 + 2.0 * phi2 + 2.0 * phi3 + phi4);

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

            n = n + 1;
        }
        if (i == 1) {
            previous_r = r;
            previous_rho = rho;
            previous_phi = phi;
            previous_p = p;
            previous_y = y;
            previous_m = m;
            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            f2 = dI;
            previousI = previousPreviousI + 0.5 * (f1 + f2) * rStep;
            I = previousI;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p1 = dp;
            m1 = dm;
            y1 = dy;
            phi1 = y;
            p = previous_p + rStepHalve * p1;
            m = previous_m + rStepHalve * m1;
            y = previous_y + rStepHalve * y1;
            phi = previous_phi + rStepHalve * phi1;

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

            r = previous_r + rStepHalve;
            r2 = r * r;
            r3 = r * r2;
            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            I = previousI + 0.5 * (dI + f2) * rStepHalve;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p2 = dp;
            m2 = dm;
            y2 = dy;
            phi2 = y;
            p = previous_p + rStepHalve * p2;
            m = previous_m + rStepHalve * m2;
            y = previous_y + rStepHalve * y2;
            phi = previous_phi + rStepHalve * phi2;

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            I = previousI + 0.5 * (dI + f2) * rStepHalve;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p3 = dp;
            m3 = dm;
            y3 = dy;
            phi3 = y;
            p = previous_p + rStep * p3;
            m = previous_m + rStep * m3;
            y = previous_y + rStep * y3;
            phi = previous_phi + rStep * phi3;

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

            r = previous_r + rStep;
            r2 = r * r;
            r3 = r * r2;
            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            I = previousI + 0.5 * (dI + f2) * rStep;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p4 = dp;
            m4 = dm;
            y4 = dy;
            phi4 = y;
            p = previous_p + rStep6 * (p1 + 2.0 * p2 + 2.0 * p3 + p4);
            m = previous_m + rStep6 * (m1 + 2.0 * m2 + 2.0 * m3 + m4);
            y = previous_y + rStep6 * (y1 + 2.0 * y2 + 2.0 * y3 + y4);
            phi = previous_phi + rStep6 * (phi1 + 2.0 * phi2 + 2.0 * phi3 + phi4);

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

            n = n + 1;
        }
    }
    while ((p / pStart) > toleranse) {
        previous_r = r;
        previous_rho = rho;
        previous_phi = phi;
        previous_p = p;
        previous_y = y;
        previous_m = m;
        yAnnen = y * y;
        phiAnnen = phi * phi;
        phiTredje = phi * phiAnnen;
        phiFjerde = phiAnnen * phiAnnen;
        dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
        I = previousPreviousI + rStep * (dI + 4.0 * f2 + f1) / 3.0;
        previousPreviousI = previousI;
        previousI = I;
        f1 = f2;
        f2 = dI;
        dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
        dm = 4.0 * pi * r2 * rho;
        dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
        p1 = dp;
        m1 = dm;
        y1 = dy;
        phi1 = y;
        p = previous_p + rStepHalve * p1;
        m = previous_m + rStepHalve * m1;
        y = previous_y + rStepHalve * y1;
        phi = previous_phi + rStepHalve * phi1;

        rhoMax = previous_rho;
        rhoHalve = rhoMax / 2.0;
        rhoMin = 0.0;
        x = xKonstant * pow(rhoMax, 1.0 / 3.0);
        pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
        pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        pMin = 0.0;
        while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
            if (pHalve > p) {
                rhoMax = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
            else {
                rhoMin = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
        }
        rho = (rhoMax + rhoMin) / 2.0;

        r = previous_r + rStepHalve;
        r2 = r * r;
        r3 = r * r2;
        yAnnen = y * y;
        phiAnnen = phi * phi;
        phiTredje = phi * phiAnnen;
        phiFjerde = phiAnnen * phiAnnen;
        dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
        I = previousI + 0.5 * (dI + f2) * rStepHalve;
        dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
        dm = 4.0 * pi * r2 * rho;
        dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
        p2 = dp;
        m2 = dm;
        y2 = dy;
        phi2 = y;
        p = previous_p + rStepHalve * p2;
        m = previous_m + rStepHalve * m2;
        y = previous_y + rStepHalve * y2;
        phi = previous_phi + rStepHalve * phi2;

        rhoMax = previous_rho;
        rhoHalve = rhoMax / 2.0;
        rhoMin = 0.0;
        x = xKonstant * pow(rhoMax, 1.0 / 3.0);
        pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
        pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        pMin = 0.0;
        while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
            if (pHalve > p) {
                rhoMax = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
            else {
                rhoMin = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
        }
        rho = (rhoMax + rhoMin) / 2.0;

        yAnnen = y * y;
        phiAnnen = phi * phi;
        phiTredje = phi * phiAnnen;
        phiFjerde = phiAnnen * phiAnnen;
        dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
        I = previousI + 0.5 * (dI + f2) * rStepHalve;
        dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
        dm = 4.0 * pi * r2 * rho;
        dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
        p3 = dp;
        m3 = dm;
        y3 = dy;
        phi3 = y;
        p = previous_p + rStep * p3;
        m = previous_m + rStep * m3;
        y = previous_y + rStep * y3;
        phi = previous_phi + rStep * phi3;

        rhoMax = previous_rho;
        rhoHalve = rhoMax / 2.0;
        rhoMin = 0.0;
        x = xKonstant * pow(rhoMax, 1.0 / 3.0);
        pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
        pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        pMin = 0.0;
        while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
            if (pHalve > p) {
                rhoMax = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
            else {
                rhoMin = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
        }
        rho = (rhoMax + rhoMin) / 2.0;

        r = previous_r + rStep;
        r2 = r * r;
        r3 = r * r2;
        yAnnen = y * y;
        phiAnnen = phi * phi;
        phiTredje = phi * phiAnnen;
        phiFjerde = phiAnnen * phiAnnen;
        dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
        I = previousI + 0.5 * (dI + f2) * rStep;
        dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
        dm = 4.0 * pi * r2 * rho;
        dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
        p4 = dp;
        m4 = dm;
        y4 = dy;
        phi4 = y;
        p = previous_p + rStep6 * (p1 + 2.0 * p2 + 2.0 * p3 + p4);
        m = previous_m + rStep6 * (m1 + 2.0 * m2 + 2.0 * m3 + m4);
        y = previous_y + rStep6 * (y1 + 2.0 * y2 + 2.0 * y3 + y4);
        phi = previous_phi + rStep6 * (phi1 + 2.0 * phi2 + 2.0 * phi3 + phi4);

        rhoMax = previous_rho;
        rhoHalve = rhoMax / 2.0;
        rhoMin = 0.0;
        x = xKonstant * pow(rhoMax, 1.0 / 3.0);
        pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
        pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        pMin = 0.0;
        while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
            if (pHalve > p) {
                rhoMax = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
            else {
                rhoMin = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
        }
        rho = (rhoMax + rhoMin) / 2.0;

        if (n == stepWhenPhiStrikes) {
            phi = faktor * phiVacuum;
        }
        n = n + 1;
    }
    phi = previous_phi;

}


void equilibriumValues(int n, 
                       double& p, 
                       double& m,
                       double& r, 
                       double& phi, 
                       double& dphi, 
                       double& rho, 
                       double pStart,
                       double toleranse, 
                       double mu, 
                       double lambda,
                       double M, 
                       double rStep, 
                       double phiVacuum,
                       int stepWhenPhiStrikes,
                       double faktor,
                       double* rVector, 
                       double* pVector,  
                       double* dpVector,
                       double* ddpVector, 
                       double* mVector, 
                       double* dmVector, 
                       double* rhoVector,
                       double* drhoVector, 
                       double* phiVector,
                       double* dphiVector,
                       double* ddphiVector,
                       double* e2betaVector,
                       double* dbetaVector,
                       double* dalphaVector,
                       double* ddalphaVector,
                       double* gVector,
                       double* fVector, 
                       double* dVvector,
                       double* ddVvector,
                       double* Vvector,
                       double* kVector, 
                       double* dkVector, 
                       double* ddrhoVector,
                       double epsilon,
                       double A,
                       double B, 
                       double xKonstant, 
                       double pKonstant) {

    double y, 
           dy, 
           dp,
           dm,
           I, 
           dI,
           mu2,
           mu4,
           M2,
           phiAnnen,
           phiTredje,
           phiFjerde,
           r2,
           r3,
           pi, 
           constant1,
           previousPreviousI, 
           previousI, 
           previous_p, 
           previous_m, 
           previous_phi,
           previous_rho, 
           previous_r,
           previous_y,
           f1,
           f2,
           y1, 
           y2, 
           y3, 
           y4,
           p1,
           p2,
           p3,
           p4,
           m1,
           m2,
           m3,
           m4,
           phi1, 
           phi2,
           phi3,
           phi4,
           rStepHalve, 
           rStep6,
           yAnnen, 
           dV,
           V,  
           f,
           g,
           term1,
           term2,
           term3,
           term4,
           term5,
           term6,
           term7,
           drho,  
           dalpha,
           ddalpha, 
           ddp, 
           rhoMax,
           rhoHalve,
           rhoMin,
           pMax,
           pHalve,
           pMin,
           x,
           x2,
           x3,
           u,
           rotU, 
           rotUUU,
           xPlussRotU, 
           v,
           k,
           dk, 
           ddrho, 
           rhoPowTredjedel;

    rStepHalve = rStep / 2.0;
    rStep6 = rStep / 6.0;
    phiAnnen = phi * phi;
    phiTredje = phi * phiAnnen;
    phiFjerde = phiAnnen * phiAnnen;
    mu2 = mu * mu;
    mu4 = mu2 * mu2;
    M2 = M * M;
    r2 = r * r;
    r3 = r * r2;
    pi = acos(-1.0L);
    constant1 = (mu2 * mu2) / (4.0 * sqrt(lambda));

    y = dphi;
    yAnnen = y * y;

    for (int i = 0; i < 2; i++) {
        if (i == 0) {
            previous_r = r;
            previous_rho = rho;
            previous_phi = phi;
            previous_p = p;
            previous_y = y;
            previous_m = m;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            f1 = dI;
            previousPreviousI = dI * r;
            I = previousPreviousI;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);

            dV = -mu2 * phi + lambda * phiTredje;
            V = -0.5 * mu2 * phiAnnen + 0.25 * lambda * phiFjerde + constant1;
            f = (1.0 - 0.5 * phiAnnen / M2) / (M2 + phiAnnen);
            g = phi / (M2 + 0.5 * phiAnnen);
            term1 = g * (3.0 * p - rho);
            term2 = 2.0 * pi * I;
            term3 = p + 0.5 * yAnnen - V;
            term4 = 4.0 * pi * r3 * term3;
            term5 = term2 + m + term4;
            term6 = r * (r - 2.0 * term2 - 2.0 * m);
            term7 = yAnnen + rho + p;

            rhoPowTredjedel = pow(rho, 1.0 / 3.0);
            x = xKonstant * rhoPowTredjedel;
            x2 = x * x;
            x3 = x2 * x;
            u = 1.0 + x2;
            v = 2.0 * x2 / 3.0 - 1.0;
            rotU = sqrt(u);
            rotUUU = sqrt(u * u * u);
            xPlussRotU = x + rotU;
            k = A * (rotU * v + x2 * v / rotU + 4.0 * x2 * rotU / 3.0 + (1.0 + x / rotU) / xPlussRotU);
            dk = A * (3.0 * x * v / rotU + 4.0 * x * rotU + 8.0 * x3 / (3.0 * rotU) - x3 * v / rotUUU + ((1.0 / rotU - x2 / rotUUU) * xPlussRotU - (1.0 + x / rotU) * (1.0 + x / rotU)) / (xPlussRotU * xPlussRotU));
            drho = dp * pow(rho, 2.0 / 3.0) / k;

            dalpha = term5 / term6;
            ddalpha = ((2.0 * pi * dI + dm + 12.0 * pi * r2 * term3 + 4.0 * pi * r3 * (y * dy - dV * y + dp)) * term6 - term5 * (term6 / r + r * (1.0 - 4.0 * pi * dI - 2.0 * dm))) / (term6 * term6);
            ddp = f * (3.0 * p - rho) * yAnnen + g * (3.0 * dp - drho) * y + term1 * dy - (2.0 * y * dy + drho + dp) * term5 / term6 - term7 * ddalpha;

            ddrho = (2.0 * drho * dp / 3.0 + rho * ddp - dk * B * drho * drho / rhoPowTredjedel) / (k * rhoPowTredjedel);

            rVector[i] = r;
            pVector[i] = p;
            mVector[i] = m;
            rhoVector[i] = rho;
            phiVector[i] = phi;
            dpVector[i] = dp;
            dmVector[i] = dm;
            drhoVector[i] = drho;
            ddrhoVector[i] = ddrho;
            dphiVector[i] = y;
            ddpVector[i] = ddp;
            ddphiVector[i] = dy;
            ddalphaVector[i] = ddalpha;
            fVector[i] = f;
            gVector[i] = g;
            Vvector[i] = V;
            dVvector[i] = dV;
            dalphaVector[i] = dalpha;
            ddVvector[i] = -mu2 + 3.0 * lambda * phiAnnen;
            e2betaVector[i] = r2 / term6;
            dbetaVector[i] = 0.5 * (1.0 / r - (1.0 - 2.0 * dm - 4.0 * pi * dI) / (r - 2.0 * m - 4.0 * pi * I));
            kVector[i] = k;
            dkVector[i] = dk;

            p1 = dp;
            m1 = dm;
            y1 = dy;
            phi1 = y;
            p = previous_p + rStepHalve * p1;
            m = previous_m + rStepHalve * m1;
            y = previous_y + rStepHalve * y1;
            phi = previous_phi + rStepHalve * phi1;

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

            r = previous_r + rStepHalve;
            r2 = r * r;
            r3 = r * r2;
            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            I = previousPreviousI + 0.5 * (dI + f1) * rStepHalve;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p2 = dp;
            m2 = dm;
            y2 = dy;
            phi2 = y;
            p = previous_p + rStepHalve * p2;
            m = previous_m + rStepHalve * m2;
            y = previous_y + rStepHalve * y2;
            phi = previous_phi + rStepHalve * phi2;

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            I = previousPreviousI + 0.5 * (dI + f1) * rStepHalve;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p3 = dp;
            m3 = dm;
            y3 = dy;
            phi3 = y;
            p = previous_p + rStep * p3;
            m = previous_m + rStep * m3;
            y = previous_y + rStep * y3;
            phi = previous_phi + rStep * phi3;

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

            r = previous_r + rStep;
            r2 = r * r;
            r3 = r * r2;
            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            I = previousPreviousI + 0.5 * (dI + f1) * rStep;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p4 = dp;
            m4 = dm;
            y4 = dy;
            phi4 = y;
            p = previous_p + rStep6 * (p1 + 2.0 * p2 + 2.0 * p3 + p4);
            m = previous_m + rStep6 * (m1 + 2.0 * m2 + 2.0 * m3 + m4);
            y = previous_y + rStep6 * (y1 + 2.0 * y2 + 2.0 * y3 + y4);
            phi = previous_phi + rStep6 * (phi1 + 2.0 * phi2 + 2.0 * phi3 + phi4);

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

        }
        if (i == 1) {
            previous_r = r;
            previous_rho = rho;
            previous_phi = phi;
            previous_p = p;
            previous_y = y;
            previous_m = m;
            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            f2 = dI;
            previousI = previousPreviousI + 0.5 * (f1 + f2) * rStep;
            I = previousI;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);

            dV = -mu2 * phi + lambda * phiTredje;
            V = -0.5 * mu2 * phiAnnen + 0.25 * lambda * phiFjerde + constant1;
            f = (1.0 - 0.5 * phiAnnen / M2) / (M2 + phiAnnen);
            g = phi / (M2 + 0.5 * phiAnnen);
            term1 = g * (3.0 * p - rho);
            term2 = 2.0 * pi * I;
            term3 = p + 0.5 * yAnnen - V;
            term4 = 4.0 * pi * r3 * term3;
            term5 = term2 + m + term4;
            term6 = r * (r - 2.0 * term2 - 2.0 * m);
            term7 = yAnnen + rho + p;

            rhoPowTredjedel = pow(rho, 1.0 / 3.0);
            x = xKonstant * rhoPowTredjedel;
            x2 = x * x;
            x3 = x2 * x;
            u = 1.0 + x2;
            v = 2.0 * x2 / 3.0 - 1.0;
            rotU = sqrt(u);
            rotUUU = sqrt(u * u * u);
            xPlussRotU = x + rotU;
            k = A * (rotU * v + x2 * v / rotU + 4.0 * x2 * rotU / 3.0 + (1.0 + x / rotU) / xPlussRotU);
            dk = A * (3.0 * x * v / rotU + 4.0 * x * rotU + 8.0 * x3 / (3.0 * rotU) - x3 * v / rotUUU + ((1.0 / rotU - x2 / rotUUU) * xPlussRotU - (1.0 + x / rotU) * (1.0 + x / rotU)) / (xPlussRotU * xPlussRotU));
            drho = dp * pow(rho, 2.0 / 3.0) / k;

            dalpha = term5 / term6;
            ddalpha = ((2.0 * pi * dI + dm + 12.0 * pi * r2 * term3 + 4.0 * pi * r3 * (y * dy - dV * y + dp)) * term6 - term5 * (term6 / r + r * (1.0 - 4.0 * pi * dI - 2.0 * dm))) / (term6 * term6);
            ddp = f * (3.0 * p - rho) * yAnnen + g * (3.0 * dp - drho) * y + term1 * dy - (2.0 * y * dy + drho + dp) * term5 / term6 - term7 * ddalpha;

            ddrho = (2.0 * drho * dp / 3.0 + rho * ddp - dk * B * drho * drho / rhoPowTredjedel) / (k * rhoPowTredjedel);

            rVector[i] = r;
            pVector[i] = p;
            mVector[i] = m;
            rhoVector[i] = rho;
            phiVector[i] = phi;
            dpVector[i] = dp;
            dmVector[i] = dm;
            drhoVector[i] = drho;
            ddrhoVector[i] = ddrho;
            dphiVector[i] = y;
            ddpVector[i] = ddp;
            ddphiVector[i] = dy;
            ddalphaVector[i] = ddalpha;
            fVector[i] = f;
            gVector[i] = g;
            Vvector[i] = V;
            dVvector[i] = dV;
            dalphaVector[i] = dalpha;
            ddVvector[i] = -mu2 + 3.0 * lambda * phiAnnen;
            e2betaVector[i] = r2 / term6;
            dbetaVector[i] = 0.5 * (1.0 / r - (1.0 - 2.0 * dm - 4.0 * pi * dI) / (r - 2.0 * m - 4.0 * pi * I));
            kVector[i] = k;
            dkVector[i] = dk;

            p1 = dp;
            m1 = dm;
            y1 = dy;
            phi1 = y;
            p = previous_p + rStepHalve * p1;
            m = previous_m + rStepHalve * m1;
            y = previous_y + rStepHalve * y1;
            phi = previous_phi + rStepHalve * phi1;

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

            r = previous_r + rStepHalve;
            r2 = r * r;
            r3 = r * r2;
            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            I = previousI + 0.5 * (dI + f2) * rStepHalve;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p2 = dp;
            m2 = dm;
            y2 = dy;
            phi2 = y;
            p = previous_p + rStepHalve * p2;
            m = previous_m + rStepHalve * m2;
            y = previous_y + rStepHalve * y2;
            phi = previous_phi + rStepHalve * phi2;

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            I = previousI + 0.5 * (dI + f2) * rStepHalve;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p3 = dp;
            m3 = dm;
            y3 = dy;
            phi3 = y;
            p = previous_p + rStep * p3;
            m = previous_m + rStep * m3;
            y = previous_y + rStep * y3;
            phi = previous_phi + rStep * phi3;

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;

            r = previous_r + rStep;
            r2 = r * r;
            r3 = r * r2;
            yAnnen = y * y;
            phiAnnen = phi * phi;
            phiTredje = phi * phiAnnen;
            phiFjerde = phiAnnen * phiAnnen;
            dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
            I = previousI + 0.5 * (dI + f2) * rStep;
            dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
            dm = 4.0 * pi * r2 * rho;
            dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
            p4 = dp;
            m4 = dm;
            y4 = dy;
            phi4 = y;
            p = previous_p + rStep6 * (p1 + 2.0 * p2 + 2.0 * p3 + p4);
            m = previous_m + rStep6 * (m1 + 2.0 * m2 + 2.0 * m3 + m4);
            y = previous_y + rStep6 * (y1 + 2.0 * y2 + 2.0 * y3 + y4);
            phi = previous_phi + rStep6 * (phi1 + 2.0 * phi2 + 2.0 * phi3 + phi4);

            rhoMax = previous_rho;
            rhoHalve = rhoMax / 2.0;
            rhoMin = 0.0;
            x = xKonstant * pow(rhoMax, 1.0 / 3.0);
            pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
            pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            pMin = 0.0;
            while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
                if (pHalve > p) {
                    rhoMax = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                    pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
                else {
                    rhoMin = rhoHalve;
                    rhoHalve = (rhoMax + rhoMin) / 2.0;
                    x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                    pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                    x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                    pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                }
            }
            rho = (rhoMax + rhoMin) / 2.0;
        }
    }
    for (int i = 2; i < n; i++) {
        previous_r = r;
        previous_rho = rho;
        previous_phi = phi;
        previous_p = p;
        previous_y = y;
        previous_m = m;
        yAnnen = y * y;
        phiAnnen = phi * phi;
        phiTredje = phi * phiAnnen;
        phiFjerde = phiAnnen * phiAnnen;
        dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
        I = previousPreviousI + rStep * (dI + 4.0 * f2 + f1) / 3.0;
        previousPreviousI = previousI;
        previousI = I;
        f1 = f2;
        f2 = dI;
        dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
        dm = 4.0 * pi * r2 * rho;
        dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);

        dV = -mu2 * phi + lambda * phiTredje;
        V = -0.5 * mu2 * phiAnnen + 0.25 * lambda * phiFjerde + constant1;
        f = (1.0 - 0.5 * phiAnnen / M2) / (M2 + phiAnnen);
        g = phi / (M2 + 0.5 * phiAnnen);
        term1 = g * (3.0 * p - rho);
        term2 = 2.0 * pi * I;
        term3 = p + 0.5 * yAnnen - V;
        term4 = 4.0 * pi * r3 * term3;
        term5 = term2 + m + term4;
        term6 = r * (r - 2.0 * term2 - 2.0 * m);
        term7 = yAnnen + rho + p;

        rhoPowTredjedel = pow(rho, 1.0 / 3.0);
        x = xKonstant * rhoPowTredjedel;
        x2 = x * x;
        x3 = x2 * x;
        u = 1.0 + x2;
        v = 2.0 * x2 / 3.0 - 1.0;
        rotU = sqrt(u);
        rotUUU = sqrt(u * u * u);
        xPlussRotU = x + rotU;
        k = A * (rotU * v + x2 * v / rotU + 4.0 * x2 * rotU / 3.0 + (1.0 + x / rotU) / xPlussRotU);
        dk = A * (3.0 * x * v / rotU + 4.0 * x * rotU + 8.0 * x3 / (3.0 * rotU) - x3 * v / rotUUU + ((1.0 / rotU - x2 / rotUUU) * xPlussRotU - (1.0 + x / rotU) * (1.0 + x / rotU)) / (xPlussRotU * xPlussRotU));
        drho = dp * pow(rho, 2.0 / 3.0) / k;

        dalpha = term5 / term6;
        ddalpha = ((2.0 * pi * dI + dm + 12.0 * pi * r2 * term3 + 4.0 * pi * r3 * (y * dy - dV * y + dp)) * term6 - term5 * (term6 / r + r * (1.0 - 4.0 * pi * dI - 2.0 * dm))) / (term6 * term6);
        ddp = f * (3.0 * p - rho) * yAnnen + g * (3.0 * dp - drho) * y + term1 * dy - (2.0 * y * dy + drho + dp) * term5 / term6 - term7 * ddalpha;

        ddrho = (2.0 * drho * dp / 3.0 + rho * ddp - dk * B * drho * drho / rhoPowTredjedel) / (k * rhoPowTredjedel);

        //ofile << setw(15) << setprecision(13) << drho << "\t";
        //ofile << setw(15) << setprecision(13) << ddrho << "\t";
        //ofile << setw(15) << setprecision(13) << r << "\t";
        //ofile << setw(15) << setprecision(13) << phi << endl;

        rVector[i] = r;
        pVector[i] = p;
        mVector[i] = m;
        rhoVector[i] = rho;
        phiVector[i] = phi;
        dpVector[i] = dp;
        dmVector[i] = dm;
        drhoVector[i] = drho;
        ddrhoVector[i] = ddrho;
        dphiVector[i] = y;
        ddpVector[i] = ddp;
        ddphiVector[i] = dy;
        ddalphaVector[i] = ddalpha;
        fVector[i] = f;
        gVector[i] = g;
        Vvector[i] = V;
        dVvector[i] = dV;
        dalphaVector[i] = dalpha;
        ddVvector[i] = -mu2 + 3.0 * lambda * phiAnnen;
        e2betaVector[i] = r2 / term6;
        dbetaVector[i] = 0.5 * (1.0 / r - (1.0 - 2.0 * dm - 4.0 * pi * dI) / (r - 2.0 * m - 4.0 * pi * I));
        kVector[i] = k;
        dkVector[i] = dk;

        p1 = dp;
        m1 = dm;
        y1 = dy;
        phi1 = y;
        p = previous_p + rStepHalve * p1;
        m = previous_m + rStepHalve * m1;
        y = previous_y + rStepHalve * y1;
        phi = previous_phi + rStepHalve * phi1;

        rhoMax = previous_rho;
        rhoHalve = rhoMax / 2.0;
        rhoMin = 0.0;
        x = xKonstant * pow(rhoMax, 1.0 / 3.0);
        pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
        pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        pMin = 0.0;
        while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
            if (pHalve > p) {
                rhoMax = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
            else {
                rhoMin = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
        }
        rho = (rhoMax + rhoMin) / 2.0;

        r = previous_r + rStepHalve;
        r2 = r * r;
        r3 = r * r2;
        yAnnen = y * y;
        phiAnnen = phi * phi;
        phiTredje = phi * phiAnnen;
        phiFjerde = phiAnnen * phiAnnen;
        dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
        I = previousI + 0.5 * (dI + f2) * rStepHalve;
        dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
        dm = 4.0 * pi * r2 * rho;
        dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
        p2 = dp;
        m2 = dm;
        y2 = dy;
        phi2 = y;
        p = previous_p + rStepHalve * p2;
        m = previous_m + rStepHalve * m2;
        y = previous_y + rStepHalve * y2;
        phi = previous_phi + rStepHalve * phi2;

        rhoMax = previous_rho;
        rhoHalve = rhoMax / 2.0;
        rhoMin = 0.0;
        x = xKonstant * pow(rhoMax, 1.0 / 3.0);
        pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
        pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        pMin = 0.0;
        while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
            if (pHalve > p) {
                rhoMax = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
            else {
                rhoMin = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
        }
        rho = (rhoMax + rhoMin) / 2.0;

        yAnnen = y * y;
        phiAnnen = phi * phi;
        phiTredje = phi * phiAnnen;
        phiFjerde = phiAnnen * phiAnnen;
        dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
        I = previousI + 0.5 * (dI + f2) * rStepHalve;
        dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
        dm = 4.0 * pi * r2 * rho;
        dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
        p3 = dp;
        m3 = dm;
        y3 = dy;
        phi3 = y;
        p = previous_p + rStep * p3;
        m = previous_m + rStep * m3;
        y = previous_y + rStep * y3;
        phi = previous_phi + rStep * phi3;

        rhoMax = previous_rho;
        rhoHalve = rhoMax / 2.0;
        rhoMin = 0.0;
        x = xKonstant * pow(rhoMax, 1.0 / 3.0);
        pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
        pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        pMin = 0.0;
        while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
            if (pHalve > p) {
                rhoMax = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
            else {
                rhoMin = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
        }
        rho = (rhoMax + rhoMin) / 2.0;

        r = previous_r + rStep;
        r2 = r * r;
        r3 = r * r2;
        yAnnen = y * y;
        phiAnnen = phi * phi;
        phiTredje = phi * phiAnnen;
        phiFjerde = phiAnnen * phiAnnen;
        dI = r2 * (yAnnen - mu2 * phiAnnen + lambda * phiFjerde / 2.0 + 2.0 * constant1);
        I = previousI + 0.5 * (dI + f2) * rStep;
        dp = (phi * (3.0 * p - rho) * y) / (M2 + phiAnnen / 2.0) - (yAnnen + rho + p) * (2 * pi * I + m + 4 * pi * r3 * (p + yAnnen / 2.0 + mu2 * phiAnnen / 2.0 - lambda * phiFjerde / 4.0 - constant1)) / (r * (r - 4 * pi * I - 2.0 * m));
        dm = 4.0 * pi * r2 * rho;
        dy = -2.0 * y / r - mu2 * phi + lambda * phiTredje - (phi * (3 * p - rho)) / (M2 + phiAnnen / 2.0);
        p4 = dp;
        m4 = dm;
        y4 = dy;
        phi4 = y;
        p = previous_p + rStep6 * (p1 + 2.0 * p2 + 2.0 * p3 + p4);
        m = previous_m + rStep6 * (m1 + 2.0 * m2 + 2.0 * m3 + m4);
        y = previous_y + rStep6 * (y1 + 2.0 * y2 + 2.0 * y3 + y4);
        phi = previous_phi + rStep6 * (phi1 + 2.0 * phi2 + 2.0 * phi3 + phi4);

        rhoMax = previous_rho;
        rhoHalve = rhoMax / 2.0;
        rhoMin = 0.0;
        x = xKonstant * pow(rhoMax, 1.0 / 3.0);
        pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
        pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
        pMin = 0.0;
        while (fabs(1.0 - rhoMin / rhoMax) > epsilon) {
            if (pHalve > p) {
                rhoMax = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMax, 1.0 / 3.0);
                pMax = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
            else {
                rhoMin = rhoHalve;
                rhoHalve = (rhoMax + rhoMin) / 2.0;
                x = xKonstant * pow(rhoMin, 1.0 / 3.0);
                pMin = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
                x = xKonstant * pow(rhoHalve, 1.0 / 3.0);
                pHalve = pKonstant * (x * sqrt(1 + x * x) * (2.0 * x * x / 3.0 - 1) + log(x + sqrt(1 + x * x)));
            }
        }
        rho = (rhoMax + rhoMin) / 2.0;

        if (i == stepWhenPhiStrikes) {
            phi = faktor * phiVacuum;
        }
    }
}

void findingAlpha(int n, 
                  double* dphiVector, 
                  double* ddphiVector,
                  double* pVector, 
                  double* dpVector,
                  double* rhoVector, 
                  double* dVvector,
                  double* rVector, 
                  double* Vvector,
                  double* mVector, 
                  double* alphaVector) {

    double presentAlpha, 
           Alpha, 
           pi, 
           R, 
           M, 
           V, 
           dp,
           dphi, 
           f1,
           f2,
           f3 , 
           ddphi, 
           dV,
           r, 
           rho,
           p, 
           previous_p,
           previousPrevious_p, 
           previous_alpha,
           previousPrevious_alpha;

    pi = acos(-1.0L);
    R = rVector[n - 1];
    M = mVector[n - 1];
    V = Vvector[n - 1];

    Alpha = 0.5 * log(1.0 - 2.0 * M / R - 8.0 * pi * R * R * V / 3.0);

    dphi = dphiVector[n - 1];
    ddphi = ddphiVector[n - 1];
    dp = dpVector[n - 1];
    p = pVector[n - 1];
    rho = rhoVector[n - 1];
    dV = dVvector[n - 1];
    r = rVector[n - 1];

    f1 = (1.0 + dphi * (ddphi - dV) / dp + 2.0 * dphi * dphi / (dp * r)) / (dphi * dphi + rho + p);

    previousPrevious_p = p;

    dphi = dphiVector[n - 2];
    ddphi = ddphiVector[n - 2];
    dp = dpVector[n - 2];
    p = pVector[n - 2];
    rho = rhoVector[n - 2];
    dV = dVvector[n - 2];
    r = rVector[n - 2];

    previous_p = p;

    f2 = (1.0 + dphi * (ddphi - dV) / dp + 2.0 * dphi * dphi / (dp * r)) / (dphi * dphi + rho + p);

    dp = pVector[n - 2] - pVector[n - 1];
    alphaVector[n - 1] = Alpha;
    alphaVector[n - 2] = Alpha - dp * 0.5 * (f1 + f2);

    previousPrevious_alpha = alphaVector[n - 1];
    previous_alpha = alphaVector[n - 2];

    for (int i = n; i > 2; i--) {
        dphi = dphiVector[i - 3];
        ddphi = ddphiVector[i - 3];
        dp = dpVector[i - 3];
        p = pVector[i - 3];
        rho = rhoVector[i - 3];
        dV = dVvector[i - 3];
        r = rVector[i - 3];

        f3 = (1.0 + dphi * (ddphi - dV) / dp + 2.0 * dphi * dphi / (dp * r)) / (dphi * dphi + rho + p);
        dp = (p - previousPrevious_p) / 2.0;

        presentAlpha = previousPrevious_alpha - dp * (f3 + 4.0 * f2 + f1) / 3.0;
        alphaVector[i - 3] = presentAlpha;
        f1 = f2;
        f2 = f3;

        previousPrevious_p = previous_p;
        previousPrevious_alpha = previous_alpha;
        previous_p = p;
        previous_alpha = presentAlpha;

    }

}

void functions(int n, 
               double* rVector,
               double* pVector,
               double* dpVector,
               double* ddpVector,
               double* mVector,
               double* dmVector,
               double* rhoVector,
               double* drhoVector, 
               double* phiVector,
               double* dphiVector, 
               double* ddphiVector,
               double* e2betaVector,
               double* dbetaVector,
               double* alphaVector,
               double* dalphaVector,
               double* ddalphaVector,
               double* gVector, 
               double* fVector,
               double* dVvector,
               double* ddVvector, 
               double* Vvector,
               double* a,
               double* c,
               double* d,
               double* aTilde,
               double* cTilde,
               double* dTilde,
               double* omegaTerm,
               double* bTerm,
               double* bBarTerm,
               double* nevner,
               double* ddrhoVector,
               double* kVector, 
               double* dkVector, 
               double B) {

    double pi, 
           r, 
           p,
           dp,
           ddp,
           m,
           dm,
           rho,
           drho,
           ddrho, 
           phi, 
           dphi,
           ddphi,
           e2beta,
           dbeta,
           alpha,
           dalpha,
           ddalpha,
           g,
           f,
           V,
           dV,
           ddV,
           dphi2,
           term1,
           term2,
           term3,
           term4,
           term5, 
           term6,
           term7,
           term8,
           term9, 
           term10,
           term11,
           term12,
           term13, 
           term14,
           k,
           dk,
           r2;

    pi = acos(-1.0L);


    for (int i = 0; i < n; i++) {
        r = rVector[i];
        p = pVector[i];
        dp = dpVector[i];
        ddp = ddpVector[i];
        m = mVector[i];
        dm = dmVector[i];
        rho = rhoVector[i];
        drho = drhoVector[i];
        ddrho = ddrhoVector[i];
        phi = phiVector[i];
        dphi = dphiVector[i];
        ddphi = ddphiVector[i];
        e2beta = e2betaVector[i];
        dbeta = dbetaVector[i];
        alpha = alphaVector[i];
        dalpha = dalphaVector[i];
        ddalpha = ddalphaVector[i];
        g = gVector[i];
        f = fVector[i];
        V = Vvector[i];
        dV = dVvector[i];
        ddV = ddVvector[i];
        k = kVector[i];
        dk = dkVector[i];

        term1 = pow(rho, 2.0 / 3.0);
        term2 = 3.0 * k / term1 - 1.0;
        term3 = 3.0 * p - rho;
        term4 = p + rho;
        term5 = 2.0 * term4 / r;
        term6 = dp + drho + term5;
        term7 = k / term1;
        term8 = 1.0 + term7;
        term9 = dbeta + dalpha;
        term10 = g * term3;
        term11 = 2.0 * dalpha + 1.0 / r;
        term12 = term7 * term4;
        term13 = dp + drho;
        term14 = term1 * term1;

        r2 = r * r;
        dphi2 = dphi * dphi;


        a[i] = -2.0 / r;
        bTerm[i] = ddV - f * term3 + g * term10 * term2;
        c[i] = g * term2 * term4;
        d[i] = g * term2 * term6;

        aTilde[i] = (-term7 * term6 - dalpha * term8 * term4 - term9 * term7 * term4 + g * dphi * term2 * term4 - dk * drho * B * term4 / term14 + 2.0 * k * drho * term4 / (3.0 * rho * term1) - term7 * term13) / term12;

        bBarTerm[i] = -dk * B * drho * term6 / term14 + 2.0 * k * drho * term6 / (3.0 * rho * term1) - term7 * (ddp + ddrho - term5 / r + 2.0 * term13 / r) - dalpha * term8 * term6 - term9 * (term7 * term6 + term4 * term11) + g * dphi * term2 * term6;

        cTilde[i] = (-term7 * term10 + 2.0 * dphi * dalpha + term9 * dphi - term10) / term12;

        dTilde[i] = (-dk * B * drho * term10 / term14 + 2.0 * k * drho * term10 / (3.0 * rho * term1) - term7 * f * dphi * term3 - term7 * g * (3.0 * dp - drho) - dalpha * term8 * term10 + term9 * (dphi * term11 - dV - term7 * term10) - f * dphi * term3 + g * term10 * dphi * term2) / term12;

        omegaTerm[i] = -e2beta * term4 / exp(2.0 * alpha);
        nevner[i] = term12;

    }

}

void eigenvalues(double stepWhenPhiStrikes, 
                 double rStep, 
                 double xiMax,
                 int n, 
                 int N,
                 double* a,
                 double* b, 
                 double* c, 
                 double* d,
                 double* aTilde,
                 double* bTilde,
                 double* cTilde,
                 double* dTilde, 
                 double* omegaTerm, 
                 double* bTerm,
                 double* bBarTerm,
                 double* xi, 
                 double* deltaPhi,
                 double* xiZero,
                 double* omega2,
                 double* nevner) {

    double h, 
           h2, 
           h3, 
           h2Halve;

    h = rStep;
    h2 = h * h;
    h3 = h * h2;
    h2Halve = h2 / 2.0;

    for (int j = 0; j < N; j++) {

        for (int i = 0; i < n; i++) {
            b[i] = bTerm[i] - omega2[j];
            bTilde[i] = (bBarTerm[i] + omega2[j] * omegaTerm[i]) / nevner[i];
        }

        xi[n - 1] = xiMax;
        xi[n - 2] = (1.0 + rStep * bTilde[n - 1] / aTilde[n - 1]) * xiMax;
        deltaPhi[n - 1] = 0.0;
        deltaPhi[n - 2] = h2Halve * (d[n - 1] - c[n - 1] * bTilde[n - 1] / aTilde[n - 1]) * xiMax;

        for (int i = n - 1; i > 1; i--) {

            xi[i - 2] = (2.0 * ((2.0 + h * a[i - 1]) * (2.0 + h2 * bTilde[i - 1]) - h3 * d[i - 1] * cTilde[i - 1]) * xi[i - 1] + ((h * aTilde[i - 1] - 2.0) * (2.0 + h * a[i - 1]) - h2 * cTilde[i - 1] * c[i - 1]) * xi[i] + 2.0 * (h2 * dTilde[i - 1] * (2.0 + h * a[i - 1]) - (2.0 * h + h3 * b[i - 1]) * cTilde[i - 1]) * deltaPhi[i - 1] + 4.0 * h * cTilde[i - 1] * deltaPhi[i]) / ((2.0 + h * aTilde[i - 1]) * (2.0 + h * a[i - 1]) - h2 * c[i - 1] * cTilde[i - 1]);

            deltaPhi[i - 2] = (2.0 * ((2.0 + h * aTilde[i - 1]) * (2.0 + h2 * b[i - 1]) - h3 * dTilde[i - 1] * c[i - 1]) * deltaPhi[i - 1] + ((h * a[i - 1] - 2.0) * (2.0 + h * aTilde[i - 1]) - h2 * cTilde[i - 1] * c[i - 1]) * deltaPhi[i] + 2.0 * (h2 * d[i - 1] * (2.0 + h * aTilde[i - 1]) - (2.0 * h + h3 * bTilde[i - 1]) * c[i - 1]) * xi[i - 1] + 4.0 * h * c[i - 1] * xi[i]) / ((2.0 + h * aTilde[i - 1]) * (2.0 + h * a[i - 1]) - h2 * c[i - 1] * cTilde[i - 1]);

            if (i <= stepWhenPhiStrikes) {
                xi[i - 2] = (2.0 * ((2.0 + h * a[i - 1]) * (2.0 + h2 * bTilde[i - 1]) - h3 * d[i - 1] * cTilde[i - 1]) * xi[i - 1] + ((h * aTilde[i - 1] - 2.0) * (2.0 + h * a[i - 1]) - h2 * cTilde[i - 1] * c[i - 1]) * xi[i]) / ((2.0 + h * aTilde[i - 1]) * (2.0 + h * a[i - 1]) - h2 * c[i - 1] * cTilde[i - 1]);
            }

        }

        xiZero[j] = xi[0];
        //ofile << setw(15) << setprecision(20) << omega2[j] << "\t";
        //ofile << setw(15) << setprecision(20) << xiZero[j] << endl;
    }
}

void tellerEgenverdier(int N,
                       int& antallEgenverdier, 
                       double* xiZero) {

    double fraction, 
           forrigeDelta, 
           Delta;

    fraction = 1.0;
    forrigeDelta = xiZero[0];
    antallEgenverdier = 0;

    // Teller antall egenverdier
    for (int i = 1; i < N; i++) {
        Delta = xiZero[i];
        fraction = Delta / forrigeDelta;
        if (fraction < 0) {
            antallEgenverdier += 1;
        }
        forrigeDelta = Delta;
    }

    //cout << "Antall egenverdier er: " << antallEgenverdier << endl;

}

void finnerEgenverdier(double stepWhenPhiStrikes, 
                       double xiMax, 
                       double rStep, 
                       int antallEgenverdier, 
                       int N, 
                       int n, 
                       double* xiZero,
                       double* omega2,
                       double* a,
                       double* b, 
                       double* c,
                       double* d,
                       double* bTerm,
                       double* aTilde, 
                       double* bTilde,
                       double* cTilde,
                       double* dTilde,
                       double* bBarTerm, 
                       double* omegaTerm,
                       double* nevner,
                       double* xi, 
                       double* deltaPhi,
                       double* egenverdiVektor,
                       double* besteXiNullVektor) {
    
    double h,  
           h2, 
           h3,
           h2Halve,
           forrigeXi,
           Xi, 
           fraction, 
           epsilon,
           besteXiNull,
           omegaMin,
           omegaMax,
           omegaStep;
    
    int teller;

    h = rStep;
    h2 = h * h;
    h3 = h * h2;
    h2Halve = h2 / 2.0;

    // Presisjon
    epsilon = 1e-13;

    // Lager matrise som skal inneholde de forskjellige OmegaMax og OmegaMin omkring hver egenverdi
    double** omegaMatrise = new double* [antallEgenverdier];

    for (int i = 0; i < antallEgenverdier; i++) {
        omegaMatrise[i] = new double[2];
    }

    // Fyller opp omegaMatrise
    forrigeXi = xiZero[0];
    teller = 0;

    for (int i = 1; i < N; i++) {
        Xi = xiZero[i];
        fraction = Xi / forrigeXi;
        if (fraction < 0) {
            omegaMatrise[teller][0] = omega2[i - 1];
            omegaMatrise[teller][1] = omega2[i];
            teller += 1;
        }
        forrigeXi = Xi;
    }

    N = 3; // Senker N for å speede opp koden. Presisjonen blir ikke dårligere siden det er epsilon som bestemmer presisjonen.
    // Lavere N gjør bare at while-løkken nedenfor entres flere ganger. Dette gjør ikke noe, siden for ved høy N, entres den allerede få ganger. Lavere N gjør at den bare entres maks et par ganger ekstra.

    // Finner egenverdiene innenfor ønsket presisjon
    for (int k = 0; k < antallEgenverdier; k++) {
        omegaMin = omegaMatrise[k][0];
        omegaMax = omegaMatrise[k][1];
        omegaStep = (omegaMax - omegaMin) / ((double)(N - 1));
        // Break for å produsere den fundamentale egenverdien. Dette sparer tid, siden vi kun er interesserte i den første egenverdien; egenverdi nummer k=0.
        if (k > 0) {
          break;
        }        
        cout << "Finner egenverdi nummer " << k << " til ønsket presisjon." <<endl;
        while (fabs(1.0 - omegaMin / omegaMax) > epsilon) {
            for (int j = 0; j < N; j++) {
                omega2[j] = omegaMin + j * omegaStep;
            }

            for (int j = 0; j < N; j++) {
                for (int i = 0; i < n; i++) {
                    b[i] = bTerm[i] - omega2[j];
                    bTilde[i] = (bBarTerm[i] + omega2[j] * omegaTerm[i]) / nevner[i];
                }

                xi[n - 1] = xiMax;
                xi[n - 2] = (1.0 + rStep * bTilde[n - 1] / aTilde[n - 1]) * xiMax;
                deltaPhi[n - 1] = 0.0;
                deltaPhi[n - 2] = h2Halve * (d[n - 1] - c[n - 1] * bTilde[n - 1] / aTilde[n - 1]) * xiMax;

                for (int i = n - 1; i > 1; i--) {

                    xi[i - 2] = (2.0 * ((2.0 + h * a[i - 1]) * (2.0 + h2 * bTilde[i - 1]) - h3 * d[i - 1] * cTilde[i - 1]) * xi[i - 1] + ((h * aTilde[i - 1] - 2.0) * (2.0 + h * a[i - 1]) - h2 * cTilde[i - 1] * c[i - 1]) * xi[i] + 2.0 * (h2 * dTilde[i - 1] * (2.0 + h * a[i - 1]) - (2.0 * h + h3 * b[i - 1]) * cTilde[i - 1]) * deltaPhi[i - 1] + 4.0 * h * cTilde[i - 1] * deltaPhi[i]) / ((2.0 + h * aTilde[i - 1]) * (2.0 + h * a[i - 1]) - h2 * c[i - 1] * cTilde[i - 1]);

                    deltaPhi[i - 2] = (2.0 * ((2.0 + h * aTilde[i - 1]) * (2.0 + h2 * b[i - 1]) - h3 * dTilde[i - 1] * c[i - 1]) * deltaPhi[i - 1] + ((h * a[i - 1] - 2.0) * (2.0 + h * aTilde[i - 1]) - h2 * cTilde[i - 1] * c[i - 1]) * deltaPhi[i] + 2.0 * (h2 * d[i - 1] * (2.0 + h * aTilde[i - 1]) - (2.0 * h + h3 * bTilde[i - 1]) * c[i - 1]) * xi[i - 1] + 4.0 * h * c[i - 1] * xi[i]) / ((2.0 + h * aTilde[i - 1]) * (2.0 + h * a[i - 1]) - h2 * c[i - 1] * cTilde[i - 1]);

                    if (i <= stepWhenPhiStrikes) {
                        xi[i - 2] = (2.0 * ((2.0 + h * a[i - 1]) * (2.0 + h2 * bTilde[i - 1]) - h3 * d[i - 1] * cTilde[i - 1]) * xi[i - 1] + ((h * aTilde[i - 1] - 2.0) * (2.0 + h * a[i - 1]) - h2 * cTilde[i - 1] * c[i - 1]) * xi[i]) / ((2.0 + h * aTilde[i - 1]) * (2.0 + h * a[i - 1]) - h2 * c[i - 1] * cTilde[i - 1]);
                    }
                }
                xiZero[j] = xi[0];
            }

            forrigeXi = xiZero[0];
            for (int j = 1; j < N; j++) {
                Xi = xiZero[j];
                fraction = Xi / forrigeXi;
                if (fraction < 0) {
                    besteXiNull = (forrigeXi + Xi) / 2.0;
                    omegaMin = omega2[j - 1];
                    omegaMax = omega2[j];
                    omegaStep = (omegaMax - omegaMin) / ((double)(N - 1));
                }
                forrigeXi = Xi;
            }
        }
        egenverdiVektor[k] = (omegaMin + omegaMax) / 2.0;
        besteXiNullVektor[k] = besteXiNull;
   
        //cout << "omegaNull2 er: " << egenverdiVektor[k] << endl;
        //ofile << setw(15) << setprecision(13) << egenverdiVektor[k] << "\t";
        //ofile << setw(15) << setprecision(13) << besteXiNull << endl;
    }

}
