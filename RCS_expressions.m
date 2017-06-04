%%Radiofrequency Circuits and Systems useful expressions
%%Version: 1.1
function RCS_expressions(option)
    remain= 1;
    choice= 0;
    while remain, 
    if nargin() > 1, 
        disp("Invalid number of arguments\n"),
    elseif nargin == 1 || choice == 1,
        if choice == 1;
            answer= answer2;
            choice= 0;
        elseif choice == 2;
            choice= 0;
        else
            answer= option;
        endif
    else
        if choice == 2;
            choice= 0;
        endif
    clc, %clean shell
    disp("________________________________________________________"), 
    disp("Radiofrequency Circuits and Systems - Useful expressions"), 
    disp("________________________________________________________"), 
    disp(" "), 
    disp("01.-  Input impedance and propagation constant of an "), 
    disp("      infinite-length transmission line"), 
    disp("02.-  T-parameter matrix of a transmission line"), 
    disp("03.-  Lossless transmission line"), 
    disp("04.-  Lossless transmission line in sinusoidal steady state"), 
    disp("05.-  Reflection coefficient"), 
    disp("06.-  Evolution of Ro along the line"), 
    disp("07.-  Voltage wave incident on the load in the general problem"), 
    disp("08.-  Voltage at a distance l respect to the load"), 
    disp("09.-  Current (towards the load) at a distance l respect"), 
    disp("      to the load"), 
    disp("10.-  Voltage standing wave ratio"), 
    disp("11.-  Lossly transmission line"), 
    disp("12.-  Attenuation"), 
    disp("13.-  Power transferred to the load"), 
    disp("14.-  Power delivered by the generator"), 
    disp("15.-  Rectangular wave guide"), 
    disp("16.-  Power density at coordinates with respect to the"), 
    disp("      transmitter"), 
    disp("17.-  Effective area"), 
    disp("18.-  Friis transmission equation"), 
    disp("19.-  Array of two antennas with phase shift"), 
    disp("20.-  Gain of a lamda/2 dipole (arm length = lamda/4)"), 
    disp(" "),
    disp("________________________________________________________"), 
    disp("Another useful expressions"),  
    disp("________________________________________________________"),
    disp("21.-  Trigonometrical arctangent"), 
    disp("22.-  L.T. length in C.C."), 
    disp("23.-  L.T. length in O.C."), 
    disp("24.-  _ViL_ = V2⁻ * delay"),
    disp(" "),
    disp("________________________________________________________"), 
    disp("Conversion magnitude"),  
    disp("________________________________________________________"),
    disp("30.-  alpha = Attenuation / 20log(e¹)"), 
    disp("31.-  mW a dBm ---> 10log(mW)"),
    disp("32.-  dBm a mW ---> 10^(dBm/10)\n"),

    answer= input("What's your menu selection? "); 
    endif;

    switch answer
    case 1, %Input Z and propagation constant of an infinite-length T.L.
        disp("EQUATION:    Z0= sqrt(Z/Y),     gamma = sqrt(Z*Y)"), 
        Z= input("Z (or L)? "); 
        Y= input("Y (or C)? "); 
        Z0= sqrt(Z/Y); 
        G= sqrt(Z*Y);
        printf("\tZ0 = %e\n\tgamma= %e\n", Z0, G);
    case 2, %T-parameter matrix of a transmission line
        disp("\n      |   s11= cosh(theta)        s12= Z0*sinh(theta)   |"),
        disp(  "  T = |                                                 |"),
        disp(  "      |   s21= 1/sinh(theta)      s22= cosh(theta)      |\n"),
    case 3, %Lossless transmission line
        disp("EQUATION:    theta= gamma*l = Tau*s = (l/Vp)*s,   Vp = 1/sqrt(L*C)"),
        unknown= input("\tCalculate theta (1) or Vp (2)? "); 
        switch unknown
        case 1,
            l= input("l? ");
            G= input("gamma? ");
            Tau= input("Tau? ");
            s= input("s? ");
            Vp= input("Vp? ");
            if Vp == 0, 
                disp("ERROR: Vp can't be zero\n"),
            else
                O1= G*l;
                printf("\ttheta = gamma * l = %e\n", O1);
                disp(""),
                O2= Tau*s;
                printf("\ttheta = Tau * s = %e\n", O2);
                disp(""),
                O3= (l/Vp)*s;
                printf("\ttheta = (l/Vp)*s = %e\n", O3);
                disp(""),
            endif
        case 2,
            Z= input("Z (or L)? "); 
            Y= input("Y (or C)? "); 
            if sqrt(Z*Y) == 0
                disp("CAUTION:    Zero division\n");
            else
                Vp= 1/(sqrt(Z*Y));
                printf("\tVp = %e\n", Vp);
            endif
        endswitch            
    case 4, %Lossless transmission line in sinusoidal steady state
        disp("EQUATION:    lamda= Vp/f"), 
        Vp= input("Vp? "); 
        f= input("f? "); 
        if f == 0, 
            disp("ERROR: f can't be zero\n"),
        else
            lamda= Vp/f;
            printf("\tlamda = %e\n", lamda);
        endif
    case 5, %Reflection coefficient
        disp("EQUATION  Ro= (Z - Z0)/(Z + Z0),       Z = Z0 * (1 + Ro)/(1 - Ro)"),
        yesno= input("\tCalculate Ro (1) or Z (2)? "); 
        switch yesno
        case 1,
            Z= input("Z? "); 
            Z0= input("Z0? "); 
            if (Z+Z0) == 0,
                disp("CAUTION:    Zero division\n");
            else
                Ro= (Z-Z0)/(Z+Z0);
                printf("\tRo = %e\n", Ro);
            endif
        case 2,
            Z0= input("Z0? "); 
            Ro= input("Ro? "); 
            if (1-Ro) == 0,
                disp("CAUTION:    Zero division\n");
            else
                Z= (1+Ro)/(1-Ro);
                printf("\tZ = %e\n", Z);
            endif
        otherwise
            error ("invalid value");
        endswitch
    case 6, %Evolution of Ro along the line
        disp("EQUATION:    Ro= RoL*exp(-2*theta)"), 
        RoL= input("RoL? ");
        alpha= input("alpha? ");
        l= input("l? ");
        theta= alpha*l;
        Ro= RoL*exp(-2*alpha*l);
        printf("INFO VALUE: theta = alpha * l = %e\n", theta);
        printf("\tRo = %e\n", Ro);
    case 7, %Voltage wave incident on the load in the general problem
        disp("EQUATION:    ViL= (1/2)*(1-RoG)*(exp(-alpha*l)/(1-RoG*RoL*exp(-2*theta)))*VG"), 
        RoG= input("RoG? ");
        RoL= input("RoL? ");
        alpha= input("alpha? ");
        l= input("l? ");
        theta= alpha*l;
        if (1-RoG*RoL*exp(-2*alpha*l)) == 0,
            disp("CAUTION:    Zero division\n");
        else
            VG= input("VG? ");
            ViL= (1/2)*(1-RoG)*(exp(-alpha*l)/(1-RoG*RoL*exp(-2*alpha*l)))*VG;
            printf("INFO VALUE: theta = alpha * l = %e\n", theta);
            printf("\tViL = %e\n", ViL);
        endif
    case 8, %Voltage at a distance l respect to the load
        disp("EQUATION:    V= ViL*(exp(theta)+RoL*exp(-theta))"), 
        RoL= input("RoL? ");
        alpha= input("alpha? ");
        l= input("l? ");
        theta= alpha*l;
        ViL= input("ViL? ");
        V= ViL*(exp(theta)+RoL*exp(-theta));
        printf("INFO VALUE: theta = alpha * l = %e\n", theta);
        printf("\tV = %e\n", V);
    case 9, %Current (towards the load) at a distance l respect to the load
        disp("EQUATION:    I= (ViL/Z0)*(exp(theta)-RoL*exp(-theta))"), 
        ViL= input("ViL? ");
        Z0= input("Z0? "); 
        if Z0 == 0, 
            disp("ERROR: Z0 can't be zero\n"),
        else
            RoL= input("RoL? ");
            alpha= input("alpha? ");
            l= input("l? ");
            theta= alpha*l;
            I= (ViL/Z0)*(exp(theta)-RoL*exp(-theta));
            printf("INFO VALUE: theta = alpha * l = %e\n", theta);
            printf("\tI = %e\n", I);
        endif
    case 10, %Voltage standing wave ratio
        disp("EQUATION:    ROE= (1+abs(Ro))/(1-abs(Ro))"), 
        Ro= input("Ro? ");
        if (1-abs(Ro)) == 0,
            disp("CAUTION:    Zero division\n");
        else
            ROE= (1+abs(Ro))/(1-abs(Ro));
            printf("\tVSWR = ROE = %e\n", ROE);
        endif
    case 11, %Lossly transmission line
        disp("EQUATION:    theta= gamma*l = alpha*l + Tau*s"),
        G= input("gamma? ");
        l= input("l? ");
        alpha= input("alpha? ");
        O1= alpha*l;
        O2= G*l;
        Tau= input("Tau? ");
        s= input("s? ");
        O3= alpha*l+Tau*s;
        printf("INFO VALUE: alpha * l = %e\n", O1);
        disp(""),
        printf("\ttheta = gamma * l = %e\n", O2);
        disp(""),
        printf("\ttheta = alpha * l + Tau * s = %e\n", O3);
        disp(""),
    case 12, %Attenuation
        disp("EQUATION:    A= alpha*20*log10(exp(1))"), 
        unknown= input("\tCalculate A (0) or alpha (1)? ");
        if unknown == 0, %General equation
            alpha= input("alpha? ");
            A= alpha*20*log10(exp(1));
            printf("\tAttenuation (dB/m) = %e\n", A);
        else
            %disp("EQUATION:    alpha = Attenuation / 20log(e¹)"), 
            A= input("A? ");
            alpha= A/(20*log10(exp(1)));
            printf("\talpha = %e\n", alpha);
        endif
    case 13, %Power transferred to the load
        disp("EQUATION:    PL= (ViL^2/(2*Z0))*(1-RoL^2)"), 
        ViL= input("ViL? ");
        Z0= input("Z0? ");
        if Z0 == 0,
            disp("ERROR: Z0 can't be zero\n"),
        else
            RoL= input("RoL? ");
            PL= (ViL^2/(2*Z0))*(1-RoL^2);
            printf("\tPL = %e [W]\n", PL);
        endif
    case 14, %Power delivered by the generator 
        disp("EQUATION:    PG= (ViL^2/(2*Z0))*(exp(2*theta)-RoL^2*exp(-2*theta))"), 
        ViL= input("ViL? ");
        Z0= input("Z0? ");
        if Z0 == 0,
            disp("ERROR: Z0 can't be zero\n"),
        else
            RoL= input("RoL? ");
            alpha= input("alpha? ");
            l= input("l? ");
            theta= alpha*l;
            printf("INFO VALUE: theta = alpha * l = %e\n", theta);
            PG= (ViL^2/(2*Z0))*(exp(2*alpha*l)-RoL^2*exp(-2*alpha*l));
            printf("\tPG = %e\n", PG);
        endif
    case 15, %Rectangular wave guide
        disp("EQUATION:    a= lamda/(2*cos(theta))"), 
        lamda= input("lamda? ");
        alpha= input("alpha? ");
        l= input("l? ");
        theta= alpha*l;
        printf("INFO VALUE: theta = alpha * l = %e\n", theta);
        if cos(theta) == 0
            disp("CAUTION:    Zero divison\n");
        else
            a= lamda/(2*cos(theta));
            printf("\t\ta = %e\n", a);
        endif
    case 16, %Power density at coordinates with respect to the transmitter
        disp("EQUATION:    P= PT/(4*pi*r^2)*GT(theta, phi)"), 
        disp("          ... revision pending ..."), 
        PT= input("PT? ");
        r= input("r? ");
        if r == 0
            disp("ERROR: r can't be zero\n");
        else
            GT= input("GT? ");
            alpha= input("alpha? ");
            l= input("l? ");
            theta= alpha*l;
            printf("INFO VALUE: theta = alpha * l = %e\n", theta);
            P= PT/(4*pi*r^2)*GT;
            printf("\tP = %e\n", P);
        endif
    case 17, %Effective area
        disp("EQUATION:    A= (lamda^2/(4*pi))*G"), 
        lamda= input("lamda? ");
        G= input("G? ");
        A= (lamda^2/(4*pi))*G;
        printf("\tA = %e\n", A);
    case 18, %Friis transmission equation
        printf("EQUATION:    PR = PT*GT*GR*(lamda/(4*pi*r))^2\n"),
        unknown= input("\tPR (0), PT (1), GT (2), GR (3), r (4), lamda (5)? ");
        if unknown == 0, %General equation
            r= input("r? ");
            if r == 0
                disp("ERROR: r can't be zero\n");
            else
                PT= input("PT? ");
                GT= input("GT? ");
                GR= input("GR? ");
                lamda= input("lamda? ");
                PR= PT*GT*GR*(lamda/(4*pi*r))^2;
                printf("\tPR = %e (W)\n", PR);
                printf("\tPR = %e (mW)\n", PR*1000);
                printf("Another expression (checking result):\n"), 
                PR= (PT*GT*GR*lamda^2)/(4*pi*r)^2;
                printf("\tPR = %e (W)\n", PR);
                printf("\tPR = %e (mW)\n", PR*1000);
            endif
        elseif unknown == 4, %Friis r isolate by two different ways
            lamda= input("lamda? ");
            PR= input("PR? "); 
            PT= input("PT? "); 
            GT= input("GT? "); 
            GR= input("GR? "); 
            r= lamda/(4*pi*sqrt(PR/(PT*GT*GR)));
            printf("\tr = %e\n", r);
            r= sqrt((PT*GT*GR*lamda^2)/((4*pi)^2*PR));
            printf("\tr = %e\n", r);
        elseif unknown == 5, %Friis lamda isolate by two different ways
            PT= input("PT? "); 
            GT= input("GT? "); 
            GR= input("GR? "); 
            if PT == 0 || GT == 0 || GR == 0,
                disp("ERROR: PT, GT, GR can't be zero\n");
            else
                PR= input("PR? "); 
                r= input("r? ");
                lamda= 4*pi*r*(sqrt(PR/(PT*GT*GR)));
                printf("\tlamda = %e\n", lamda);
                lamda= sqrt((PR*(4*pi*r)^2)/(PT*GT*GR));
                printf("\tlamda = %e\n", lamda);
            endif
        else %Friis PT, GT, GR isolate by two different ways
            lamda= input("lamda? ");
            r= input("r? ");
            if unknown == 1,
                GT= input("GT? "); 
                GR= input("GR? "); 
                PxGx= GT*GR;
                PG= "PT";
            elseif unknown == 2, 
                PT= input("PT? "); 
                GR= input("GR? ");
                PxGx= PT*GR;
                PG= "GT";
            else %unknown == 3, 
                PT= input("PT? "); 
                GT= input("GT? ");
                PxGx= PT*GT;
                PG= "GR";
            endif
            if lamda == 0 || r == 0 || GT == 0 || GR == 0, 
                disp("ERROR: lamda, r, GT, GR can't be zero\n");
            else
                PR= input("PR? ");

                GP= PR/(PxGx*(lamda/(4*pi*r))^2);
                printf("\t%s = %e\n", PG, GP);
                printf("Another expression (checking result):\n"), 
                GP= (((4*pi*r)^2)*PR)/(PxGx*lamda^2);
                printf("\t%s = %e\n", PG, GP);
            endif
        endif
    case 19, %Array of two antennas with phase shift
        disp("EQUATION:    FA= 1+exp(j*alpha)*exp(j*k*d*cos(theta))"), 
        alpha= input("alpha? ");
        k= input("k? ");
        d= input("d? ");
        theta= input("theta? ");
        printf("INFO VALUE: theta = alpha * l = %e\n", theta);
        FA= 1+exp(j*alpha)*exp(j*k*d*cos(theta));
        printf("\tFA = %e\n", FA);
    case 20, %Gain of a lamda/2 dipole (arm length = lamda/4)
        disp("EQUATION:    G= 1.64*(cos((pi/2)*cos(theta))/sin(theta))^2"), 
        alpha= input("alpha? ");
        l= input("l? ");
        theta= alpha*l;
        printf("INFO VALUE: theta = alpha * l = %e\n", theta);
        G= 1.64*(cos((pi/2)*cos(theta))/sin(theta))^2;
        printf("\tG = %e\n", G);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%    Another useful expressions      %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 21, 
        disp("EQUATION:    l= (atan(Z0/Zin)*Vp)/(2*pi*f)"), 
        Vp= input("Vp? "); 
        f= input("f? "); 
        Z0= input("Z0? "); 
        Zin= input("Zin? "); 
        l= (atan(Z0/Zin)*Vp)/(2*pi*f);
        printf("\tl = %e\n", l);
    case 22, %Curt Circuit
        Z= input("Z (or L)? "); 
        Y= input("Y (or C)? "); 
        f= input("f? "); 
        Zin= input("Zin? "); 
        Vp= 1/sqrt(Z*Y);
        printf("\tVp = %e\n", Vp);
        Z0= j*sqrt(Z/Y);
        printf("\tZ0 = %e\n", Z0);
        Tw= (2*pi*f)/Vp;
        printf("\tTw = %e\n", Tw);
        l= (atan(Zin/Z0)/Tw);
        printf("\tl = %e\n", l);
        l2= l+(Vp/f)/2 %if l < 0
        printf("\tl2 = %e\n", l2);
    case 23, %Circuit Obert
        Z= input("Z (or L)? "); 
        Y= input("Y (or C)? "); 
        f= input("f? "); 
        Zin= input("Zin? "); 
        Vp= 1/sqrt(Z*Y);
        printf("\tVp = %e\n", Vp);
        Z0= -j*sqrt(Z/Y);
        printf("\tZ0 = %e\n", Z0);
        Tw= (2*pi*f)/Vp;
        printf("\tTw = %e\n", Tw);
        l= (atan(Z0/Zin)/Tw);
        printf("\tl = %e\n", l);
    case 24,
        disp("EQUATION:    _ViL_ = V2n*exp(-theta)"),
        V2n= input("V2⁻? "); 
        alpha= input("alpha? ");
        l= input("l? ");
        theta= alpha*l;
        ViL_= V2n*exp(-theta);
        printf("\t_ViL_ = %e\n", ViL_);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%    Conversion magnitude      %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 31, %mW to dBm
        disp("mW to dBm ---> 10log(mW)"), 
        mW= input("mw? ");
        dBm= 10*log10(mW);
        printf("\tdBm= %e\n", dBm);
        remain = 0;
    case 32, %dBm to mW
        disp("dBm to mW ---> 10^(dBm/10)"), 
        dBm= input("dBm? ");
        mW= 10^(dBm/10);
        printf("\tmW= %e\n", mW);

    otherwise
        error ("invalid value");
    endswitch

    disp("________________________________________________________"),
    choice= input("\n\tNumber's menu (1)\tMenu (2)\tExit (3)\nChoose the option: "); 
    switch choice, 
    case 1, 
        answer2 = input("What's your number menu selection? "); 
        remain = 1;
    case 2,
        choice = 2;
        remain = 1;
    case 3,
        remain = 0;
        disp("\t\n\t*****   See you again. Bye!!!   *****\n"),

%{
        yesno= input("Do you want to reset variables (y/n)? ", "s"); %"s" reference to string
        switch yesno, 
        case {"Yes" "yes" "YES" "y" "Y"}, %"s" reference to string
            clear,
            remain = 1,
        case {"No" "no" "NO" "n" "N"}, 
            remain = 1;
        otherwise
            error ("invalid value");
        endswitch, 
%}

    otherwise
      error ("invalid value");
    endswitch

    endwhile
    
endfunction


%{
NOTATION TYPES
real: %f    scientific: %e  integer: %i     string: %s      
ACCURACY (examples)
real: %.1f  scientific: %.2e

ANNOTATIONS
The disp function displays the value of a variable (scalar, vector, matrix, 
string, etc.) in the same way as simply typing the name of the variable does 
and you can select the output format too.
> format long
> disp(log(10))
%}
