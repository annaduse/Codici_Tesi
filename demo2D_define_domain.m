
function demo2D_define_domain(example)

%--------------------------------------------------------------------------
% OBJECT
%--------------------------------------------------------------------------
% This demo shows:
% 1. How to define the domains.
% 2. How to compute a cubature rule with a certain algebraic degree of
%    exactness.
% 3. How to plot domain and nodes.
%--------------------------------------------------------------------------
% INPUT
%--------------------------------------------------------------------------
% example: number that selects a domain and an example.
%--------------------------------------------------------------------------
% DATE
%--------------------------------------------------------------------------
% June 15, 2023 (A. Sommariva)
%--------------------------------------------------------------------------
% COPYRIGHT
%--------------------------------------------------------------------------
% Copyright (C) 2023-
%
% Authors:
% Alvise Sommariva, Marco Vianello.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------


if nargin < 1, example=1; end

switch example

    case 1

        domain_str='polygon';
        M=129;
        th=linspace(0,2*pi,M);
        v=[cos(th').*(sin(2*th')) sin(th').*(sin(2*th'))];
        vertices=v(1:end-1,:);
        domain_structure=define_domain(domain_str,'vertices',vertices);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 2

        domain_str='disk';
        center=[0 0];
        radius=1;
        domain_structure=define_domain(domain_str,'center',center,...
            'radius',radius);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 3

        domain_str='lune';
        centers=[0 0; 1 0];
        radii=[0.8 0.7];
        domain_structure=define_domain(domain_str,'centers',centers,...
            'radii',radii);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);


    case 4

        domain_str='circular-annular-sector';
        center=[1 1.5];
        radii=[0.5 2];
        angles=[pi/4; -pi/2];
        domain_structure=define_domain(domain_str,'center',center,...
            'radii',radii,'angles',angles);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);


    case 5

        domain_str='sector';

        center=[1 1.5];
        radius=0.5;
        angles=[pi/4; -pi/2];

        domain_structure=define_domain(domain_str,'center',center,...
            'radius',radius,'angles',angles);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 6

        domain_str='asymmetric-circular-sector';

        centers=[1 1; 2 2];
        radius=3;
        angles=[pi/4; -pi/2];

        domain_structure=define_domain(domain_str,'centers',centers,...
            'radius',radius,'angles',angles);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);
        nodes=XW(:,1:2);

        plot_2D(domain_structure,[],nodes);

    case 7

        domain_str='asymmetric-annulus';

        centers=[1 1; 2 2];
        radii=[3 1];

        domain_structure=define_domain(domain_str,'centers',centers,...
            'radii',radii);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 8

        domain_str='vertical-circular-zone';

        center=[1 1];
        radius=1;
        angles=[pi/6 -pi/2+pi/6];

        domain_structure=define_domain(domain_str,'center',center,...
            'radius',radius);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 9

        domain_str='horizontal-circular-zone';

        center=[0 0];
        radius=1;
        angles=[-pi/4 -pi/6];

        domain_structure=define_domain(domain_str,'center',center,...
            'radius',radius);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 10

        domain_str='circular-segment';

        center=[1 1];
        radius=2;
        angles=[-pi/4 -pi/2];

        domain_structure=define_domain(domain_str,'center',center,...
            'radius',radius,'angles',angles);


        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);


    case 11

        domain_str='symmetric-lens';

        center=1;
        radius=2.5;

        domain_structure=define_domain(domain_str,'center',center,...
            'radius',radius);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 12

        domain_str='butterfly';

        center=[1 2];
        radius=1;
        angles=[-pi/4 pi/6]; % angles(1) < angles(2)

        domain_structure=define_domain(domain_str,'center',center,...
            'radius',radius,'angles',angles);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 13 % ISSUES

        domain_str='candy';

        center=0.5;
        radius=1;
        angle=-1.5;

        domain_structure=define_domain(domain_str,'center',center,...
            'radius',radius,'angle',angle);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 14

        domain_str='NURBS';

        P=[1 0; 1 1; 0 1;  -1 0;  -1 -1; -0.4 0.6; 0 0.5; 1 -1;  1 0];
        order=4;
        L=size(P,1); tt=linspace(0,1,L-order+2);
        ttinit=zeros(1,order); ttmid=tt(2:end-1);  ttfin=ones(1,order);
        knots=[ttinit ttmid ttfin];
        c=1/sqrt(2); w=[1 c/2 1 c/4 1 c/4 1 c 1];

        geometry_NURBS=makeNURBSarc('free','P',...
            P,'knots',knots,'weights',w,'order',order);

        domain_structure=define_domain(domain_str,...
            'geometry',geometry_NURBS);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 15

        domain_str='union-disks';

        disks=[4.172e-01 6.4437e-01 9.630e-01
            4.965e-02 3.786e-01 5.468e-01
            9.027e-01 8.115e-01 5.211e-01];

        centers=disks(:,1:2);
        radii=disks(:,3);

        domain_structure=define_domain(domain_str,'centers',centers,...
            'radii',radii);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 16

        domain_str='square';

        domain_structure=define_domain(domain_str);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 17

        domain_str='unit-square[0,1]x[0,1]';

        domain_structure=define_domain(domain_str);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 18

        domain_str='rectangle';

        xLimits=[0.5 1];
        yLimits=[2 3];

        domain_structure=define_domain(domain_str,'xLimits',xLimits,...
            'yLimits',yLimits);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 19

        domain_str='triangle';

        vertices=[0 0; 1 0; 1 1; 0 0];

        domain_structure=define_domain(domain_str,'vertices',vertices);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 20

        domain_str='polygcirc';

        extrema=[0.25 0; 0 0.25];
        vertices=[0.5 0; 0.5 0.5; 0 0.5];
        center=[0 0];
        radius=0.25;
        convex=0;

        domain_structure=define_domain(domain_str,'extrema',extrema,...
            'vertices',vertices,'center',center,'radius',radius,...
            'convex',convex);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);

    case 21

        domain_str='unit-simplex';

        domain_structure=define_domain(domain_str);

        deg=10; flag_compression=1;
        XW=define_cub_rule(domain_structure,deg,flag_compression);

        nodes=XW(:,1:2);
        plot_2D(domain_structure,[],nodes);


end

