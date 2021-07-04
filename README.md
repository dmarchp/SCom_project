# Complex Systems Project
Simulation of a SIR epidemic model on any kind of network. The code is divided in different modules: SIR simulation routines, network input reading routines and network properties routines.

## Reading and setting up the network:
The module **set_up_net** is found in read_net.f90. The subroutine *read_net* reads the network specified in an input file. This input file must contain all the links of the network in columns, with no duplication (e.g. 1 2, 1 3, 2 3) or duplicated (e.g. 1 2, 1 3, 2 1, 2 3, 3 1, 3 2). Global module variables will be initialized when reading the net, which can be used for posteior applications:
- *Integers*:
   - num_links: number of links in the system
   - Nnodes: number of nodes of the network
   - min_node, max_node: index of the min/max node of the input network
   - min_degree, max_degree: min/max degree among all nodes
- *Allocated vectors while reading the net*:
   - degree(min_node:max_node): contains the degree of each node of the network
   - neighbors(2\*num_links): contains the neighbors of each node, on increasing node index
   - p_ini, p_fin (min_node:max_node): initial and final position of the neigbors of node *i* in the vector neighbors
 
## Computing network properties:
The module **net_props** is found in net_props.f90 (makes use of *set_up_net*). Here can be found many subroutines to compute the input network properties: degree distribution, average nearest neighbor degree distribution, number of triangles and the degree dependent or clustering spectrum.
 
## SIR network dynamics:
The module **sir** is found in sir.f90 (makes use of *set_up_net*). Here are defined global variables and routines necessary for the SIR evolution. The important global variables are:
- *Integers*:
   - parameters: SUSC = 0, INF = 1, RECOV = 2 (parameters to designate node states)
   - Nactive_links, Ninfected: number of active or infectious links in the network (i.e. from an infected to a susceptible) and number of infected nodes.
- *Allocatable Arrays*:
   - node_state(min_node:max_node): state of every node in the network
   - infected_nodes(min_node:max_node): list current infected nodes
   - active_links(2,num_links): list of all current active links. (1,:)-> infected nodes, (2,:)-> susceptible nodes
   - neighbors_SIR(2\*num_links): contains the position of each link defined in *neighbors* in the array *active_links*. 0 if the link is not active.
    
Before the dynamics starts, the file input_SIR.txt is read. Here the infection/recovering rates can be specified, as long as the initial population fraction on each of the SIR states. The infection rate is *lambda* and the recovering rate is *delta*. It is recomended to set up the time scale of the dynamics to one of these variables, for instance leaving *delta=1* and varying only its counterpart *lambda*.
The SIR dynamics function as Gillespie type algorithm. The main routine that contains the evolution loop is *SIR_evolution*. It calls every other routine, in the following order:
 - *compute_reac_probs*: From the input reaction rates and the number of infected/active links computes the reaction probability of each reaction to happen in    the following step (infection or recovery).
 - *compute_reac_time*: From the reaction probabilities returns the time until next reaction, assuming it is exponentialy distributed.
 - *choose_reac*: From the reaction probabilities it choses the reaction that will occur at the just computed time, from a simple tower sampling.
 - *update_system*: According to the chosen reaction it will update the system state.
    - Infection: one of the active links is chosen randomly, the susceptible node is infected. Then the link is removed from *active_links* (as all other active links that pointed to the newly infected node), the just infected node is added to *infected_nodes* and new active links are searched.
     - Recovery: on of the infected nodes is chosen randomly, is is removed from *infected_nodes* and its state is set to RECOV, so it can not be susceptible anymore. Active links that emerged from this node are removed.
        
        
    
    

