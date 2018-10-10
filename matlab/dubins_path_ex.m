
function [] = dubins_path_ex()

    global LSL LSR RSL RSR RLR LRL TSL TSR LST RST TST LSEG SSEG RSEG TSEG INF TYPENAME DIRNAME DIRDATA;

    dubins_const;

%     global LSL LSR RSL RSR RLR LRL TSL TSR LST RST TST LSEG SSEG RSEG TSEG INF;
%     
%     LSL = 0;
%     LSR = 1;
%     RSL = 2;
%     RSR = 3;
%     RLR = 4;
%     LRL = 5;
%     
%     TSL = 6;
%     TSR = 7;
%     LST = 8;
%     RST = 9;
%     TST = 10;
%     
%     LSEG = 0;
%     SSEG = 1;
%     RSEG = 2;
%     TSEG = 3;
%     
%     INF=100000.0;
%     
%     TYPENAME = ['LSL'; 'LSR'; 'RSL'; 'RSR'; 'RLR'; 'LRL'; 'TSL'; 'TSR'; 'LST'; 'RST'; 'TST'];
%     DIRNAME = ['LSEG'; 'SSEG'; 'RSEG'; 'TSEG'];
%     
%     DIRDATA = [LSEG, SSEG, LSEG; LSEG, SSEG, RSEG; RSEG, SSEG, LSEG; RSEG, SSEG, RSEG; RSEG, LSEG, RSEG; LSEG, RSEG, LSEG; ...
%         TSEG, SSEG, LSEG; TSEG, SSEG, RSEG; LSEG, SSEG, TSEG; RSEG, SSEG, TSEG; TSEG, SSEG, TSEG;];

    t = [INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF];
    p = [INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF];
    q = [INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF];
    cost = [INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF];
    gamma1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    gamma3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    sgn1 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
    sgn3 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        
    iRes = -1; %10;
    
    %a = [6.0, 2.0, 2.722];
    %a = [1.0, 1.0, 1.0];
    %a = [550.0, 102.0, 0.0];
    a = [546.0, 114.0, 2.024600];
    
    %b = [6.0, 0.0, 3.142];
    %b = [16.0, 2.0, 2.0];
    %b = [546.0, 114.0, 2.024600];
    b = [550.0, 102.0, 0.0];
    fprintf(1, 'a: x=%f y=%f alpha=%f\n', a(1), a(2), a(3)); 
    fprintf(1, 'b: x=%f y=%f beta=%f\n', b(1), b(2), b(3)); 
    
    a(3) = wrapTo2Pi(a(3));
    b(3) = wrapTo2Pi(b(3));
    
    vel = 15.0; %4.0;
    omega_turn = 15.0;
    omega_arc = 0.2; %1.6;
	rho = vel/omega_arc;

    H=0.05;
    %H=0.052360;
    fprintf(1, 'rho=%f H=%f\n', rho, H); 
    
    dx = b(1) - a(1);
    dy = b(2) - a(2);
    D = sqrt(dx*dx + dy*dy);
    d = D/rho;
    fprintf(1, 'dx=%f dy=%f D=%f d=%f\n', dx, dy, D, d);    
    
    theta = wrapTo2Pi(atan2(dy, dx));
    alpha = wrapTo2Pi(a(3) - theta);
    beta  = wrapTo2Pi(b(3) - theta);
    fprintf(1, 'theta=%f alpha=%f beta=%f\n', theta, alpha, beta);
    
    sa = sin(alpha);            
    sb = sin(beta);             
    ca = cos(alpha);           
    cb = cos(beta);             
    cab = cos(alpha - beta);

    %LSL
    i = LSL + 1;
    tmp0 = d+sa-sb;
    tmp2 = 2 + (d*d) - (2*cab) + (2*d*(sa - sb));
    if( tmp2 >= 0 && tmp0 >= 0) 
        tmp1 = atan((cb-ca)/tmp0);
        t(i) = wrapTo2Pi(-alpha + tmp1);
        p(i) = sqrt(tmp2);
        q(i) = wrapTo2Pi(beta - tmp1);
        %cost(i) = (t(i) + p(i) + q(i))/vel;
        cost(i) = get_cost(t(i), p(i), q(i), vel);
    end
    fprintf(1, 'LSL (%d): t=%f p=%f q=%f cost=%f\n', LSL, t(i), p(i), q(i), cost(i));
    
    %LSR
    i = LSR + 1;
    tmp1 = -2 + (d*d) + (2*cab) + (2*d*(sa+sb));
    if( tmp1 >= 0 ) 
        p(i)  = sqrt(tmp1);
        tmp2 = atan((-ca-cb)/(d+sa+sb)) - atan(-2.0/p(i));
        t(i) = wrapTo2Pi(-alpha + tmp2);
        q(i) = wrapTo2Pi(-wrapTo2Pi(beta) + tmp2);
        cost(i) = get_cost(t(i), p(i), q(i), vel);
    end    
    fprintf(1, 'LSR (%d): t=%f p=%f q=%f cost=%f\n', LSR, t(i), p(i), q(i), cost(i));
    
    %RSL
    i = RSL + 1;
    tmp1 = (d*d) - 2 + (2*cab) - (2*d*(sa+sb));
    if( tmp1 > 0 ) 
        p(i) = sqrt(tmp1);
        tmp2 = atan((ca+cb)/(d-sa-sb)) - atan(2.0/p(i));
        t(i) = wrapTo2Pi(alpha - tmp2);
        q(i) = wrapTo2Pi(beta - tmp2);
        cost(i) = get_cost(t(i), p(i), q(i), vel);
    end
    fprintf(1, 'RSL (%d): t=%f p=%f q=%f cost=%f\n', RSL, t(i), p(i), q(i), cost(i));
        
    %RSR
    i = RSR + 1;
    tmp0 = d-sa+sb;
    tmp2 = 2 + (d*d) - (2*cab) + (2*d*(sb-sa));
    if( tmp2 >= 0 && tmp0 >= 0) 
        tmp1 = atan( (ca-cb)/tmp0 );
        t(i) = wrapTo2Pi(alpha - tmp1);
        p(i) = sqrt(tmp2);
        q(i) = wrapTo2Pi(-beta + tmp1);
        cost(i) = get_cost(t(i), p(i), q(i), vel);
    end
    fprintf(1, 'RSR (%d): t=%f p=%f q=%f cost=%f\n', RSR, t(i), p(i), q(i), cost(i));
    
    %RLR
    i = RLR + 1;
    tmp1 = (6.0 - d*d + 2*cab + 2*d*(sa-sb))/8.0;
    if( abs(tmp1) < 1) 
        p(i) = acos(tmp1);
        t(i) = wrapTo2Pi(alpha - atan2(ca-cb, d-sa+sb) + wrapTo2Pi(p(i)/2.0));
        q(i) = wrapTo2Pi(alpha - beta - t(i) + wrapTo2Pi(p(i)));
        cost(i) = get_cost(t(i), p(i), q(i), vel);
    end
    fprintf(1, 'RLR (%d): t=%f p=%f q=%f cost=%f\n', RLR, t(i), p(i), q(i), cost(i));
    
    %LRL
    i = LRL + 1;
    tmp1 = (6.0 - d*d + 2*cab + 2*d*(-sa + sb))/8.0;
    if( abs(tmp1) < 1) 
        p(i) = wrapTo2Pi(acos(tmp1));
        t(i) = wrapTo2Pi(-alpha - atan2(ca-cb, d+sa-sb) + p(i)/2.0);
        q(i) = wrapTo2Pi(wrapTo2Pi(beta) - alpha -t(i) + wrapTo2Pi(p(i)));
        cost(i) = get_cost(t(i), p(i), q(i), vel);
    end
    fprintf(1, 'LRL (%d): t=%f p=%f q=%f cost=%f\n', LRL, t(i), p(i), q(i), cost(i));
    
    %TSL
    i = TSL + 1;
    tmp1 = d*d - 2*d*sb;
    if (tmp1 >= 0)
        p(i) = sqrt(tmp1);
        gamma1(i) = wrapTo2Pi(asin(-(d - sb - p(i)*cb)/(p(i)*p(i) + 1)));
        [delta, sgn1(i)] = turn_in_place(alpha, gamma1(i));
        t(i) = delta;
        q(i) = 2*pi - gamma1(i) + beta;
        cost(i) = get_cost_ex(t(i), p(i), q(i), vel, omega_turn, TSL);
    end;
    fprintf(1, 'TSL (%d): t=%f p=%f q=%f cost=%f gamma1=%f sgn1=%d\n', TSL, t(i), p(i), q(i), cost(i), gamma1(i), sgn1(i));
   
    %TSR
    i = TSR + 1;
    tmp1 = d*d + 2*d*sb;
    if (tmp1 >= 0)
        p(i) = sqrt(tmp1);
        gamma1(i) = wrapTo2Pi(asin((d + sb - p(i)*cb)/(p(i)*p(i) + 1)));
        [delta, sgn1(i)] = turn_in_place(alpha, gamma1(i));
        t(i) = delta;
        q(i) = 2*pi + gamma1(i) - beta;
        cost(i) = get_cost_ex(t(i), p(i), q(i), vel, omega_turn, TSR);
    end;
    fprintf(1, 'TSR (%d): t=%f p=%f q=%f cost=%f gamma1=%f sgn1=%d\n', TSR, t(i), p(i), q(i), cost(i), gamma1(i), sgn1(i));
    
    %LST
    i = LST + 1;
    tmp1 = d*d + 2*d*sa;
    if (tmp1 >= 0)
        p(i) = sqrt(tmp1);
        gamma3(i) = wrapTo2Pi(asin((d + sa - p(i)*ca)/(p(i)*p(i) + 1)));
        [delta, sgn3(i)] = turn_in_place(gamma3(i), beta);
        q(i) = delta;
        t(i) = 2*pi + gamma3(i) - alpha;
        cost(i) = get_cost_ex(t(i), p(i), q(i), vel, omega_turn, LST);
    end;
    fprintf(1, 'LST (%d): t=%f p=%f q=%f cost=%f gamma3=%f sgn3=%d\n', LST, t(i), p(i), q(i), cost(i), gamma3(i), sgn3(i));
    
    %RST
    i = RST + 1;
    tmp1 = d*d - 2*d*sa;
    if (tmp1 >= 0)
        p(i) = sqrt(tmp1);
        gamma3(i) = wrapTo2Pi(asin(-(d - sa - p(i)*ca)/(p(i)*p(i) + 1)));
        [delta, sgn3(i)] = turn_in_place(gamma3(i), beta);
        q(i) = delta;
        t(i) = 2*pi - gamma3(i) + alpha;
        cost(i) = get_cost_ex(t(i), p(i), q(i), vel, omega_turn, RST);
    end
    fprintf(1, 'RST (%d): t=%f p=%f q=%f cost=%f gamma3=%f sgn3=%d\n', RST, t(i), p(i), q(i), cost(i), gamma3(i), sgn3(i));
    
    %TST
    i = TST + 1;
    p(i) = d;
    gamma1(i) = 0.0;
    gamma3(i) = 0.0;
    [t(i), sgn1(i)] = turn_in_place(alpha, 0.0);
    [q(i), sgn3(i)] = turn_in_place(0.0, beta);
    cost(i) = get_cost_ex(t(i), p(i), q(i), vel, omega_turn, TST);
    fprintf(1, 'TST (%d): t=%f p=%f q=%f cost=%f gamma1=%f gamma3=%f sgn1=%d sgn3=%d\n', TST, t(i), p(i), q(i), cost(i), gamma1(i), gamma3(i), sgn1(i), sgn3(i));
  
    % Extract the best cost path
    [costMin, iMin] = min(cost); 
    fprintf(1, 'iMin=%d costMin=%f\n', iMin, costMin);
    
    if (iRes >= 0)
        iMin=iRes + 1;  
    end
    dubins_type = iMin-1;
    fprintf(1, 'dubins_type = %d\n', dubins_type);

    costMin=cost(iMin);
    gamma1Min = gamma1(iMin);
    gamma3Min = gamma3(iMin);
    sgn1Min = sgn1(iMin);
    sgn3Min = sgn3(iMin);
    tMin = t(iMin);
    qMin = q(iMin);
    length = (t(iMin) + p(iMin) + q(iMin))*rho;
    n=floor(length/H) + 1;
    fprintf(1, '\nType=%s (%d) costMin=%f length=%f n=%d gamma1Min=%f gamma3Min=%f sgn1Min=%d sgn3Min=%d\n', ... 
            TYPENAME(iMin,:), dubins_type, costMin, length, n, gamma1Min, gamma3Min, sgn1Min, sgn3Min);
    
    if (costMin == INF) 
        fprintf(1, 'Required Dubins path not found\n');
        return;
    end;
    
    i1 = DIRDATA(iMin, 1) + 1;
    i2 = DIRDATA(iMin, 2) + 1;
    i3 = DIRDATA(iMin, 3) + 1;
    fprintf(1, 'Path segments: %s (%d), %s (%d), %s (%d)\n', DIRNAME(i1,:), DIRDATA(iMin, 1), DIRNAME(i2,:), DIRDATA(iMin, 2), DIRNAME(i3,:), DIRDATA(iMin, 3));
        
    %a0 = [0.0, 0.0, a(3)];
    a0 = [0.0, 0.0, alpha];
    %fprintf(1, 'a0: %f %f %f\n', a0(1), a0(2), a0(3));
    
    p1 = t(iMin);
    p2 = p(iMin);
    p3 = q(iMin);
    fprintf(1, 'p1=%f p2=%f p3=%f\n', p1, p2, p3);
    
    q10 = dubins_segment(p1, a0, DIRDATA(iMin, 1), sgn1Min);
    fprintf(1, 'q10: %f %f %f\n', q10(1), q10(2), q10(3));
    q1 = transform(q10, a, rho, theta);
    fprintf(1, 'q1: %f %f %f\n', q1(1), q1(2), q1(3)); 
    
    q20 = dubins_segment(p2, q10, DIRDATA(iMin, 2), sgn1Min);
    fprintf(1, 'q20: %f %f %f\n', q20(1), q20(2), q20(3));
    q2 = transform(q20, a, rho, theta);
    fprintf(1, 'q2: %f %f %f\n', q2(1), q2(2), q2(3)); 
    
    q30 = dubins_segment(p3, q20, DIRDATA(iMin, 3), sgn1Min);
    %fprintf(1, 'q30: %f %f %f\n', q30(1), q30(2), q30(3));
    q3 = transform(q30, a, rho, theta);
    fprintf(1, 'q3: %f %f %f\n', q3(1), q3(2), q3(3)); 
            
    tt= 0.0;
    i=1;
    while (tt < length)
        
       tt0 = tt/rho;
       
       if (tt0 < p1) 
           type = DIRDATA(iMin, 1);
           qt = dubins_segment(tt0, a0, type, sgn1Min);    
       elseif (tt0 < (p1+p2)) 
           type = DIRDATA(iMin, 2);
           qt = dubins_segment(tt0-p1, q10, type, 0);
       else 
           type = DIRDATA(iMin, 3);
           qt = dubins_segment(tt0-p1-p2, q20, type, sgn3Min); 
       end;       
       qt = transform(qt, a, rho, theta);
       
       %fprintf(1, 'i=%d tt=%f x=%f y=%f theta=%f type=%d\n', i, tt, qt(1), qt(2), qt(3), type); 
       x(i) = qt(1);
       y(i) = qt(2);
       i = i+1;
       tt  = tt + H;
    end;
    %fprintf(1, 'i=%d tt=%f\n', i, tt);
    
    plot_dubins(a, b, q1, q2, x, y, H, gamma1Min, gamma3Min, sgn1Min, sgn3Min, iMin-1, tMin, qMin);
        
    function [qt] = dubins_segment(t, qi, type, sign) 
       if (type == LSEG)
           qt(1) = qi(1) + sin(qi(3)+t) - sin(qi(3));
           qt(2) = qi(2) - cos(qi(3)+t) + cos(qi(3));
           qt(3) = wrapTo2Pi(qi(3) + t);
       elseif (type == RSEG)
           qt(1) = qi(1) - sin(qi(3)-t) + sin(qi(3));
           qt(2) = qi(2) + cos(qi(3)-t) - cos(qi(3));
           qt(3) = wrapTo2Pi(qi(3) - t);
       elseif (type == SSEG) 
           qt(1) = qi(1) + cos(qi(3)) * t;
           qt(2) = qi(2) + sin(qi(3)) * t;
           qt(3) = wrapTo2Pi(qi(3));
       elseif (type == TSEG) 
           %fprintf(1, 't=%f qi(3)=%f sign=%d type=%d\n', t, qi(3), sign, type);
           qt(1) = qi(1);
           qt(2) = qi(2);
           qt(3) = wrapTo2Pi(qi(3) + t*sign);
           %fprintf(1, 'qt: %f %f %f\n', qt(1), qt(2), qt(3));
       end;
    end

    function [u] = transform_xu(x, x0, rho, theta)
        dx0 = x(1) - x0(1);
        dy0 = x(2) - x0(2);
        u(1) = (dx0*cos(theta) + dy0*sin(theta))*rho;
        u(2) = (-dx0*sin(theta) + dy0*cos(theta))*rho;
        u(3) = wrapTo2Pi(x(3) - theta);
    end
    
    function [x] = transform(u, x0, rho, theta)
         x(1) = (u(1)*cos(theta) - u(2)*sin(theta))*rho + x0(1);
         x(2) = (u(1)*sin(theta) + u(2)*cos(theta))*rho + x0(2);
         x(3) = wrapTo2Pi(u(3) + theta);
    end


    function [turn_angle, turn_sgn] = turn_in_place(alpha, gamma)
     
        turn_angle = abs(gamma - alpha);
        %fprintf(1, 'turn_angle = %f\n', turn_angle);
        if (turn_angle > pi) 
            turn_angle = abs(2*pi - turn_angle);
            %fprintf(1, 'turn_angle = %f\n', turn_angle);
        end;
        
        if (0 <= alpha && alpha < pi)
            if (alpha <= gamma && gamma < alpha + pi)
                turn_sgn=1;
            else
                turn_sgn=-1;
            end;
        else
            if (alpha-pi <= gamma && gamma < alpha)
                turn_sgn=-1;
            else
                turn_sgn=1;
            end;
        end;
        %fprintf(1, 'turn_sgn = %d\n', turn_sgn);
    end

    function [cost] = get_cost(s1, s2, s3, v)
        cost = (s1 + s2 + s3)/v;
    end

    function [tm] = segment_time(s, v, omega, seg_type)
        if (seg_type == TSEG)
            tm = s/omega;
        else
            tm = s/v;
        end
    end

    function [cost] = get_cost_ex(s1, s2, s3, v, omega, path_type)
        
        cost =  segment_time(s1, v, omega, DIRDATA(path_type, 1)) + ... 
                segment_time(s2, v, omega, DIRDATA(path_type, 2)) + ... 
                segment_time(s3, v, omega, DIRDATA(path_type, 3));
    end
end