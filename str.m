function alpha = str(Aa,Ad,lambda_a,lambda_d,k)
    %coupling strength function which returns strengths for ascending or
    %descending connections
    if (k < 0)
        alpha = Aa*exp(-abs(k)/lambda_a);
    elseif (k > 0)
        alpha = Ad*exp(-abs(k)/lambda_d);
    else
        alpha = 1;  %Is either 0 or 1 depending on the model
        %alpha = 0; %1 for NeuralModel, 0 for Phase Model
    end
end