%test Space class and methods
%Parallel to B0
sample_space = Space(4, 0.05, 0);
%test BOLD calculation
sample_space.plot_bold(0.5, 0.5);

%test BOLD and dpf calculation
sample_space.plot_df(0.5, 0.5);
sample_space.plot_df(1, 0.5);


%perpendicular to B0
sample_space = Space(4, 0.05, pi / 2);
%test BOLD calculation
sample_space.plot_bold(0.5, 0.5);

%test BOLD and dpf calculation
sample_space.plot_df(0.5, 0.5);
sample_space.plot_df(1, 0.5);