#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
#include "stepperdopr5.h"

using namespace std;

/*Solves equation with
 
	theta' = phi
	phi' = -alpha*sin(theta) - epsilon((theta/theta_o)^2 - 1)*phi + eta(theta1-theta2)
 */

struct Pend {
    Doub alpha;
    Doub theta_o;
    Doub epsilon;
    Doub other_theta;
    Doub eta;
    
    Pend(Doub alpha_in, Doub theta_o_in, Doub epsilon_in, Doub eta_in, Doub other_theta_in) {
        alpha = alpha_in;
        theta_o = theta_o_in;
        epsilon = epsilon_in;
        eta = eta_in;
        other_theta = other_theta_in;
    };
    
    //calculates ODE for whichever metrenome it's focused on
    //holds other angle constant during all calculations
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx){
        const Doub theta = y[0];
        const Doub phi = y[1];
        
        Doub thetaFrc = (theta/theta_o) * (theta/theta_o);
        
        dydx[0] = phi;
        dydx[1] = -alpha * sin(theta) - epsilon * (thetaFrc - 1.0) * phi + eta*(theta - other_theta);
    }

    void update_other_theta(Doub other_in){
        other_theta = other_in;
    }
 
};

int main() {

    ofstream outfile;
    outfile.open("partE.dat");
    outfile.setf(ios::left);
    outfile << setw(16) << "# time " <<
    setw(16) << "theta1" << setw(16) << "theta2" << endl;
    outfile << "#====================================================" << endl;
    
    //vector size for holding two variables
    int nvar = 2;
    //sets absolute and relative tolerances
    //also sets minimum step size and first guess of step size
    const Doub atol = 1.e-5, rtol = atol, h1=0.001, hmin=0.0;
    
    //amount of times interval should happen
    const Doub t_in = 0.0;
    const Doub t_final = 1500.0;
    
    //random epsilon (can be altered by coder at any point)
    const Doub epsilon = 0.1;
    
    //total amount of steps
    const int N_steps = 1000;
    
    //how much t should be incremented by each iteration
    //is dependent on total amount of steps
    const Doub delta_t = (t_final - t_in)/N_steps;
    
    //vectors that holds the starting y values
    //(will be altered as code runs)
    VecDoub ystart1(nvar);
    VecDoub ystart2(nvar);
    
    //need output object to control saving of intermediate values
    Output out;

    //arbitrary constants
    const Doub alpha = 1.0;
    const Doub theta_o = 0.1;
    const Doub eta = 0.01;
    
    //initial conditions for both metronomes, can be reset by user
    ystart1[0] = 0.1;
    ystart1[1] = 0.1;

    ystart2[0] = 0.0;
    ystart2[1] = 0.0;
    

    //makes both oscillating objects
    Pend oscillator1(alpha, theta_o, epsilon, eta, ystart2[0]);
    Pend oscillator2(alpha, theta_o, epsilon, eta, ystart1[0]);


    for (Doub t = t_in; t < t_final; t+= delta_t) {
        //do Odeint, this will automatically update ystart
        //update both metrenomes
        Odeint<StepperDopr5<Pend> > ode1(ystart1,t,t+delta_t,atol,
                                            rtol,h1,hmin,out,oscillator1);
        ode1.integrate();

        Odeint<StepperDopr5<Pend> > ode2(ystart2,t,t+delta_t,atol,
                                            rtol,h1,hmin,out,oscillator2);
        ode2.integrate();

        //since the other angle is held constant during calculation,
        //needed to update manually at end of for looop
        oscillator1.update_other_theta(ystart2[0]);
        oscillator2.update_other_theta(ystart1[0]);



        outfile << setw(16) << t << setw(16) << ystart1[0];
        outfile << setw(16) << ystart2[0] << endl;


        
    }

    outfile.close();
    
    return 0;
}