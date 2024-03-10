%%
%% Demo comparing sparse and non sparse coding
% sparse tunnin is just exp(N*x) where x is the non sparse tunning and N
% the number of units

close all

angle = 0:0.01:(2*pi);

Nunits = [2 3 5 10 20]; % number of units
NN = length(Nunits);

Phases = cell(size(Nunits));           % phases of each unit

X = cell(size(Nunits));             % activation of all units    
Xs = cell(size(Nunits));       % sparse activation of all units

for inu=1:NN
    X{inu} = zeros(length(angle), Nunits(inu));
    Phases{inu} = zeros(length(Nunits(inu)),1);

    for j=1:Nunits(inu)
        if ( Nunits(inu)>2)
            % the phases of each unit are distributed evenly along the circle
            Phases{inu}(j) = pi + (j-1) * 2*pi / Nunits(inu);
        else
            % the 2 unit case is special
            % units need to be in cuadrature otherwise they are not
            % decodable
            Phases{inu}(j) = (j-1) * pi/2;
        end
        
        % the activation of each unit is just a sinusoid, shifted by the
        % phase and elevated so it is always positive. Not critical but
        % more realistic if talking about neurons
        X{inu}(:,j) = 1 + cos( angle + Phases{inu}(j) );
    end
   

    % the sparse coding is just the exponential of the normal coding
    % the exponential includes a multiplicative term that makes the tuning
    % narrower the more units are present
    Xs{inu} = exp( Nunits(inu) * X{inu} ) / exp(Nunits(inu)*2);
end

Plot(X, Phases, Nunits)
Plot(Xs, Phases, Nunits)

%%
% Change in shape of 2D manifold embedded in 3D with different tunings
Tunnings = logspace(1.5,-0.5,10);

figure
hold
for k=1:length(Tunnings)

    Tunning=Tunnings(k);

    angle = 0:0.01:(2*pi);
    X = zeros(length(angle), 3);
    for j=1:3
        Phase = (j-1)*2*pi/3;
        X(:,j)= exp(sin(angle+Phase)/Tunning);
        X(:,j) = X(:,j)/mean(X(:,j));
    end


    plot3(X(:,1),X(:,2),X(:,3),'o')
end
view(45,45)
title('manifolds of 3 ring units with different tunning curves')


function Plot(X,Phases, Nunits)

NN = length(Phases);
figure
for i= 1:NN
    subplot(5,NN,i)
    imagesc(X{i}');
    title(sprintf('%d units', Nunits(i)))
    set(gca,'xticklabel',[],'yticklabel',[])
    ylabel('units');
    xlabel('Angle')

    subplot(5,NN,i+1*NN,'nextplot','add')
    title('units and mean of units')
    plot(X{i})
    plot(mean(X{i},2),'linewidth',2,'color','k');
    set(gca,'xticklabel',[],'yticklabel',[])
    xlabel('Angle');

    subplot(5,NN,i+2*NN)
    plot(-atan2d(mean((X{i}-1).*repmat(sin(Phases{i}),height(X{i}),1),2), mean((X{i}-1).*repmat(cos(Phases{i}),height(X{i}),1),2)))  
    title('decoded angle')


    C= cov(X{i});
    M= mean(X{i});
    [V,D]= eig(C);
    P = V * diag(sqrt(1./(diag(D) + 0.1))) ;
    W1 = bsxfun(@minus, X{i}, M);
    W = W1 * P;

    subplot(5,NN,i+3*NN)
        plot(W(:,end),W(:,end-1),'o')
        set(gca,'xlim',[-2 2],'ylim',[-2 2])
        set(gca,'xticklabel',[],'yticklabel',[])
        xlabel('PCA1');
        ylabel('PCA2');

    subplot(5,NN,i+4*NN)
    if (Nunits(i)>2)
        plot(W(:,end-1),W(:,end-2),'o')
        set(gca,'xlim',[-2 2],'ylim',[-2 2])
        set(gca,'xticklabel',[],'yticklabel',[])
        xlabel('PCA2');
        ylabel('PCA3');
    end


end
end

