%% getMaxTerm.m
% Wynter Alford
% January 2022
%
% Using a rule of thumb, picks a maximum term in the Magnus expansion to
% compute up to given the length of a pulse sequence.  Outputs are based on
% what I have tended to observe in my simulations.  This set of outputs are
% designed for quickly running the code on a personal laptop; when run on a
% cluster, more terms can be computed.  This has not been configured for
% pulse sequences longer than 48 pulses in length; a 96 pulse sequence
% would most likely be much, much slower.

function maxTerm = getMaxTerm(sequenceName, plus)
    seq = getSequence(sequenceName);
    ll = length(seq.Pulses);
    if ll<5
        maxTerm = 30; % Cluster: 70
    elseif ll<11
        maxTerm = 10; % Cluster: 32
    elseif ll<25
        maxTerm = 6;  % Cluster: 9
    else
        maxTerm = 2;  % Cluster: 4
    end

    if exist('plus','var')
        maxTerm = maxTerm + plus;
    end
end