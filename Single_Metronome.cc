#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
#include "stepperdopr5.h"

using namespace std;

/*Solves equation with
 
	theta' = phi
	phi' = -alpha*sin(theta) - epsilon((theta/theta_o)^2 - 1)*phi
 */

struct Pend {
    Doub alpha;
    Doub theta_o;
    Doub epsilon;
    
    Pend(Doub alpha_in, Doub theta_o_in, Doub epsilon_in) {
        alpha = alpha_in;
        theta_o = theta_o_in;
        epsilon = epsilon_in;
    };
    
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx){
        const Doub theta = y[0];
        const Doub phi = y[1];
        
        Doub thetaFrc = (theta/theta_o) * (theta/theta_o);
        
        dydx[0] = phi;
        dydx[1] = -alpha * sin(theta) - epsilon * (thetaFrc - 1.0) * phi;
    }
    
};

int main() {

    ofstream outfile;
    outfile.open("partCAlpha1.dat");
    outfile.setf(ios::left);
    outfile << setw(16) << "# time " <<
    setw(16) << "theta" << endl;
    outfile << "#====================================================" << endl;
    
    //vector size for holding two variables
    int nvar = 2;
    //sets absolute and relative tolerances
    //also sets minimum step size and first guess of step size
    const Doub atol = 1.e-5, rtol = atol, h1=0.001, hmin=0.0;
    
    //amount of times interval should happen
    const Doub t_in = 0.0;
    const Doub t_final = 1000.0;
    
    //random epsilon (can be altered by coder at any point)
    const Doub epsilon = 0.1;
    
    //total amount of steps
    const int N_steps = 1000;
    
    //how much t should be incremented by each iteration
    //is dependent on total amount of steps
    const Doub delta_t = (t_final - t_in)/N_steps;
    
    //vector that holds the starting y values (will be altered as code runs)
    VecDoub ystart(nvar);
    
    //need output object to control saving of intermediate values
    Output out;

    //arbitrary constants
    const Doub alpha = 1.0;
    const Doub theta_o = 0.1;
    
    //initial conditions, can be reset by user
    ystart[0] = 0.0;
    ystart[1] = 0.1;
    
    Pend oscillator(alpha, theta_o, epsilon);
    for (Doub t = t_in; t < t_final; t+= delta_t) {
        Doub lastVal = ystart[0];
        //do Odeint, this will automatically update ystart
        Odeint<StepperDopr5<Pend> > ode(ystart,t,t+delta_t,atol,
                                            rtol,h1,hmin,out,oscillator);
        ode.integrate();
        outfile << setw(16) << t << setw(16) << ystart[0];

        if(lastVal < 0.0 &&  ystart[0] >= 0.0){
            //marks period;
            outfile << setw(16) << 0.0 << endl;
        }
        else{
            outfile << endl;
        }
    }

    outfile.close();
    
    
    return 0;
}