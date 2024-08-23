%% practiceDTpfe
% function to assist on practicing discrete-time partial fraction expansion
%
% inputs:
%   order: the desired order of the rational function Y(z) to be generated.
%       This can be a scalar, but can also be a vector of length 2,
%       which gives the minimum and maximum order and this is randomly
%       sampled in this range according to a uniform distribution.
%       Optional argument, default [2,4].
%   repeatedPoles: if 0: repeated poles might occur by chance, if 1:
%       repeated poles are enforced to happen.
%       Optional argument, default: 0.
%   complexPoles: if 0: there are no complex conjugate pairs of poles, if 1:
%       complex conjugate pairs of poles are enforced to happen.
%       Optional argument, default: 0.
%   num: Fix the numerator polynomial, rendering the first three
%       arguments irrelevant. The 5th argument must be provided as well.
%   den: Fix the denominator polynomial, rendering the first three
%       arguments irrelevant. The 4th argument must be provided as well.
%
% outputs:
%   PFEinfo: An mx3 matrix, where $m$ denotes the number of
%       poles of the rational function Y(z)/z. (This equals the number of
%       elementary terms in the PFE of Y(z)/z.) Each row should contain:
%       [the pole location p, its multiplicity j, its coefficient A]
%
% When running the function, various info and details are displayed:
%   - the proper rational function Y(z),
%   - the form of the PFE of Y(z)/z showing its elementary terms,
%   - the result of the PFE for Y(z)/z, as well as for Y(z),
%   - the corresponding time sequence y[k].
% In the result, poles at z=0 give rise to (possibly delayed versions of) the
% impulse sequence delta[k].
%
function [PFEinfo] = MM_Bonus(order,repeatedPoles,complexPoles,num,den)
if(nargin==5)
    % if the numerator and denominator are provided we don't have to
    % generate them. As we do not want a non-causal system, we check on this:
    if(length(num)>length(den)), error('The transfer function must be proper'); 
    else
        % find poles and zeros from given num and den 
        [z,p,k] = residue(num,den);
        res = z;
        poles = p;

    end
elseif(nargin==4)
    error('If you provide a numerator, you should also provide a denominator');
else
    % we will have to generate a numerator and denominator. We first will
    % handle the optional arguments by specifying default values.
    % (We are just providing an example start here - feel free to modify ):
   switch nargin
        case 0      % no input
            complexPoles=0;
            repeatedPoles=0;
            order=[2,4];
        case 1      % order is given 
            complexPoles=0;
            repeatedPoles=0;
        case 2      % order and repeatedPoles are given 
        complexPoles=0;
        case 3      % order, repeatedPoles and complexPoles are given
   
   end
end


if(length(order) == 2) % case when order contains min and max order values
    order = round(unifrnd(order(1),order(2)));    % uniformly choose order
end

% initialize the vectors 
zeros = ones(order,1);
poles = ones(order,1);

% randomly pick zeros and poles from unit circle of the complex plane. To
% make sure its stable
for c = 1:order 
   zeros(c) = unifrnd(-1,1);
   poles(c) = unifrnd(-1,1);
end
    
% order, repeatedPoles and complexPoles are given

% create order many poles and zeros if 'repeatedPoles' = 0 and 'complexPoles' = 0
    if(repeatedPoles == 0 && complexPoles == 0)
     
    % create order many poles and zeros if 'repeatedPoles' = 1 and 'complexPoles' = 0
    elseif(repeatedPoles == 1 && complexPoles == 0) 
        poles(2)= poles(1);
    % create order many poles and zeros if 'repeatedPoles' = 0 and 'complexPoles' = 1 
    elseif(repeatedPoles == 0 && complexPoles == 1)
         % randomly generate a complex conjugate pair of poles
         re = unifrnd(-1,1); % real part
         im = unifrnd(-1,1); % non-zero imaginary part
         poles(1:2) = [re + 1i*im, re - 1i*im]; % add the complex conjugate pair to the array of poles
    % create order many poles and zeros if 'repeatedPoles' = 1 and 'complexPoles' = 1     
    elseif(repeatedPoles == 1 && complexPoles == 1)
        if(order<4) error(['For a transfer function to have two repeated conjugate pairs in the poles of the system, ' ...
          'then the transfer function must have at least a fourth-order denominator polynomial']);
        else
         % randomly generate a complex conjugate pair of poles
         re = unifrnd(-1,1); % real part
         im = unifrnd(-1,1); % non-zero imaginary part
         poles(1:2) = [re + 1i*im, re - 1i*im]; % add the complex conjugate pair to the array of poles
         % randomply select a pole to repeat
         random_index = randi(length(poles));
         pole_to_repeat = poles(random_index);
         % case when picked a non complex pole 
         if(isreal(pole_to_repeat))
             poles(random_index+1) = pole_to_repeat;
         % case when complex pole is picked 
         else
         % copy a complex conjugate pair of poles
         re = real(pole_to_repeat); % real part
         im = imag(pole_to_repeat); % imaginary part
         % add the complex conjugate pair to the array of poles
         poles(1:4) = [re + 1i*im, re - 1i*im, re + 1i*im, re - 1i*im]; 
         end
        end
    else error(['repeatedPoles and complexPoles can only have values 1 or 0']);
    end

    num = poly(zeros);   % get the numerator polynomial from generated zeros 
    den = poly(poles);   % get the denominator polynomial from generated poles 

% PFE
% residue to find residue and k 
[res, poles, k] = residue(num,den);

% compute the multiplicity and coefficient of each pole
m = length(poles);
multiplicity = ones(m, 1);
i = 1;
while i <= m
    if imag(poles(i)) == 0 % real pole
        j = 1;
        % check for repeated pole
        while i+j <= m && abs(poles(i+j)-poles(i)) < eps 
            j = j+1;
        end
        multiplicity(i) = j;
        i = i+j;
    else % complex conjugate poles
        j = 1;
        % check for repeated pole
        while i+j <= m && abs(abs(poles(i+j))-abs(poles(i))) < eps && abs(imag(poles(i+j)) - imag(poles(i))) < eps 
            j = j+1;
        end
        multiplicity(i) = j;
        i = i+j;
    end
end

% construct PFEinfo matrix
PFEinfo = [poles(:), multiplicity(:), res(:)];

% fancy display 

disp('Proper rational function Y(z):');
fprintf('\n');
disp('Y(z) = ');
fprintf('\n');
printFraction(num,den,poles,multiplicity,0);
fprintf('\n');

% form of the PFE of Y(z)/z showing its elementary terms
disp('PFE of Y(z)/z showing its elementary terms:');
fprintf('\n');
disp('Y(z)');
disp('----- = ')
disp('  z')
fprintf('\n');
printFraction(num,den,poles,multiplicity,1);
fprintf('\n');

% the result of the PFE for Y(z)/z, as well as for Y(z)
disp('PFE for Y(z):')
fprintf('\n');
F = getPartialFractionExpansion(res, poles, k, multiplicity);
disp(F);

% Print the PFEinfo matrix
disp('PFEinfo matrix: ')
fprintf('\n');
disp('|poles|multiplicity|residue')
fprintf('\n');
disp(PFEinfo);


% invers z-transform 
syms k delta
disp('Inverse z-transform: ')
fprintf('\n');
disp('y[k] = ans')
invz = my_iztrans(F);
disp(invz);

end

%%

% Function to print the fraction 

function printFraction(numerator, denominator, poles, multiplicity, factorized)
 switch factorized
     case 0
        % Calculate the degree of the numerator and denominator
        degree_num = length(numerator) - 1;
        degree_den = length(denominator) - 1;
        
        % Print out the numerator
        fprintf('%dz^%d', numerator(1), degree_num);
        for i = 2:length(numerator)
            if numerator(i) > 0
                fprintf(' + %.4dz^%.4d', numerator(i), degree_num-i+1);
            elseif numerator(i) < 0
                fprintf(' - %.4dz^%.4d', -numerator(i), degree_num-i+1);
            end
        end
        
        % Calculate the number of dashes needed for the division line
        line_length = max(length(numerator), length(denominator))*15+2;
        
        % Print out the division line
        fprintf('\n');
        for i = 1:line_length
            fprintf('-');
        end
        fprintf('\n');
        
        % Print out the denominator
        fprintf('%dz^%d', denominator(1), degree_den);
        for i = 2:length(denominator)
            if denominator(i) > 0
                fprintf(' + %.4dz^%.4d', denominator(i), degree_den-i+1);
            elseif denominator(i) < 0
                fprintf(' - %.4dz^%.4d', -denominator(i), degree_den-i+1);
            end
        end
        fprintf('\n');

     case 1
         % Calculate the degree of the numerator and denominator
        degree_num = length(numerator) - 1;
        degree_den = length(denominator) - 1;
        
        % Print out the numerator
        fprintf('%dz^%d', numerator(1), degree_num);
        for i = 2:length(numerator)
            if numerator(i) > 0
                fprintf(' + %.4dz^%.4d', numerator(i), degree_num-i+1);
            elseif numerator(i) < 0
                fprintf(' - %.4dz^%.4d', -numerator(i), degree_num-i+1);
            end
        end
        
        % Calculate the number of dashes needed for the division line
        line_length = max(length(numerator), length(denominator))*15+2;
        
        % Print out the division line
        fprintf('\n');
        for i = 1:line_length
            fprintf('-');
        end
        fprintf('\n');
        
        % Print out the denominator
        printFactors(poles, multiplicity)
        fprintf('\n');
 end
end

%%

% Method to print out the factorized polynomial 

function printFactors(poles, multiplicity)
% Print out the factors
for i = 1:length(poles)
    if i == 1
        if multiplicity(i) == 1
            fprintf('z(z - %.4d)', poles(i));
        else 
            fprintf('z(z - %.4d)^2', poles(i));
        end
    else
        if multiplicity(i) == 1
            fprintf('*(z - %.4d)', poles(i));
        else
            fprintf('*(z - %.4d)^2', poles(i));
        end
    end
end
fprintf('\n');
end
%%

% Method to build PFE expresion 

function syms_expr = getPartialFractionExpansion(residues, poles, k, multiplicity)
    % initialize expression as an empty string
    expr = "";
    
    % loop through each term and build the expression string
    for i = 1:length(residues)
        % build residue term
        residue_term = sprintf('%.4f', residues(i));
        
        % build pole term
        if imag(poles(i)) == 0 % real pole
            if multiplicity(i) == 1
                pole_term = sprintf('/(z %+.4f)', -poles(i));
            else
                pole_term = sprintf('/(z %+.4f)^2', -poles(i));
            end
        else % complex pole
            if multiplicity(i) == 1
                pole_term = sprintf('/(z %+.4f + %.4fi)', -real(poles(i)), imag(poles(i)));
            else
                pole_term = sprintf('/(z %+.4f + %.4fi)^2', -real(poles(i)), imag(poles(i)));
            end
        end
        
        % build full term by concatenating residue and pole terms
        full_term = strcat(residue_term, pole_term);
        
        % add plus sign if there are more terms
        if i < length(residues)
            full_term = strcat(full_term, ' + ');
        end
        
        % add full term to expression string
        expr = strcat(expr, full_term);
    end
    
    % add constant term to expression string
    if k ~= 0
        const_term = sprintf(' + %.4f', k);
        expr = strcat(expr, const_term);
    end
    
    % convert expression string to a symbolic expression
    syms_expr = str2sym(expr);
end
%%

% Function to output the inverse z-transform 
function result = my_iztrans(F)
    syms k delta;
    F = vpa(iztrans(F, 'k'));
    F = subs(F, kroneckerDelta(k, 0), delta);
    result = vpa(F, 4);
end