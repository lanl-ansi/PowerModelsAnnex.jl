"""
This function inputs the ref dictionary produced by the PowerModels parser &
checks whether the network is islanded (based on the branch + dc_line "status"
entries).

The output is a new dictionary which has the same data, split per island of the 
network.

04/21/18 [cjc] 
This function is not currently maintained and may be out of date.
Use of PowerModels.connected_components (added in v0.6.0) is suitable alternative
"""
function network_islanding(ref)
  #--
  function buspair_parameters(arcs_from, branches, buses)
    buspair_indexes = collect(Set([(i,j) for (l,i,j) in arcs_from]))

    bp_angmin = Dict([(bp, -Inf) for bp in buspair_indexes])
    bp_angmax = Dict([(bp, Inf) for bp in buspair_indexes])
    bp_branch = Dict([(bp, Inf) for bp in buspair_indexes])

    for (l,branch) in branches
      i = branch["f_bus"]
      j = branch["t_bus"]

      bp_angmin[(i,j)] = max(bp_angmin[(i,j)], branch["angmin"])
      bp_angmax[(i,j)] = min(bp_angmax[(i,j)], branch["angmax"])
      bp_branch[(i,j)] = min(bp_branch[(i,j)], l)
    end

    buspairs = Dict([((i,j), Dict(
    "branch"=>bp_branch[(i,j)],
    "angmin"=>bp_angmin[(i,j)],
    "angmax"=>bp_angmax[(i,j)],
    "rate_a"=>branches[bp_branch[(i,j)]]["rate_a"],
    "tap"=>branches[bp_branch[(i,j)]]["tap"],
    "vm_fr_min"=>buses[i]["vmin"],
    "vm_fr_max"=>buses[i]["vmax"],
    "vm_to_min"=>buses[j]["vmin"],
    "vm_to_max"=>buses[j]["vmax"]
    )) for (i,j) in buspair_indexes])

    return buspairs
  end
  #--
  function biggest_generator(gens)
    biggest_gen = nothing
    biggest_value = -Inf
    for (k,gen) in gens
      if gen["pmax"] > biggest_value
        biggest_gen = gen
        biggest_value = gen["pmax"]
      end
    end
    assert(biggest_gen != nothing)
    return biggest_gen
  end
  #--------------------------------#



  #--
  #--------------------#
  been_there=zeros(Int64,1)
  to_go=zeros(Int64,1)
  mapped=zeros(Int64,1)
  my_exit_flag=0



  Isles=Dict()

  isl_cnt=0

  new_isle_flag=1
  rock_flag=0

  start_key=Int64
  start_key= first(keys(ref[:ref_buses]))
  #--------------------#




  while my_exit_flag==0


    if new_isle_flag==1

      new_isle_flag=0
      rock_flag=0

      been_there=ones(Int64,1)*start_key
      bus_arcs = ref[:bus_arcs][start_key]
      bus_arcs_dc=ref[:bus_arcs_dc][start_key]
      bus_nodes=zeros(Int64,(length(bus_arcs)+length(bus_arcs_dc)))



      #-
      cnt=0
      for (idx,jdx,kdx) in bus_arcs
        if ref[:branch][idx]["br_status"]==1
          cnt+=1
          bus_nodes[cnt]=kdx
        end
      end
      #--

      #-
      for (idx,jdx,kdx) in bus_arcs_dc
        if ref[:dc_line][idx]["br_status"]==1
          cnt+=1
          bus_nodes[cnt]=kdx
        end
      end
      #--



      if cnt>=1


        #-
        for idx in 1:cnt

          if (bus_nodes[idx] in to_go) !=true

            if idx==1
              to_go=ones(Int64,1)*bus_nodes[idx]
            else
              append!(to_go,bus_nodes[idx])
            end

          end
        end
        #--

      else

        rock_flag=1

      end

    end

    if rock_flag==0



      new_trip=shift!(to_go)
      append!(been_there,new_trip)

      bus_arcs = ref[:bus_arcs][new_trip]
      bus_arcs_dc=ref[:bus_arcs_dc][new_trip]
      bus_nodes=zeros(Int64,(length(bus_arcs)+length(bus_arcs_dc)))



      cnt=0
      for (idx,jdx,kdx) in bus_arcs
        if ref[:branch][idx]["br_status"]==1
          cnt=cnt+1
          bus_nodes[cnt]=kdx
        end
      end
      #--

      #-
      for (idx,jdx,kdx) in bus_arcs_dc
        if ref[:dcline][idx]["br_status"]==1
          cnt+=1
          bus_nodes[cnt]=kdx
        end
      end
      #--

      #-
      for idx in 1:cnt

        if (bus_nodes[idx] in to_go) !=true  && (bus_nodes[idx] in been_there) !=true
          append!(to_go,bus_nodes[idx])
        end
      end
      #--

      if length(to_go)==0


        isl_cnt+=1

        Isles[isl_cnt]=Dict()
        Isles[isl_cnt][:map]=been_there

        if isl_cnt==1
          mapped=been_there
        else
          mapped=union(mapped,been_there)
        end


        unmapped=zeros(Int64,1)
        unmapped=setdiff(keys(ref[:bus]),mapped)

        if length(unmapped)==0
          my_exit_flag=42

        else

          new_isle_flag=1
          start_key=unmapped[1]

        end
        #--


      end
      #--


    else


      isl_cnt+=1

      Isles[isl_cnt]=Dict()
      Isles[isl_cnt][:map]=been_there

      if isl_cnt==1
        mapped=been_there
      else
        mapped=union(mapped,been_there)
      end


      unmapped=zeros(Int64,1)
      unmapped=setdiff(keys(ref[:bus]),mapped)

      if length(unmapped)==0
        my_exit_flag=42

      else

        new_isle_flag=1
        start_key=unmapped[1]

      end
      #-

    end

  end



  for i in keys(Isles)



    Isles[i][:name]=ref[:name]
    Isles[i][:areas]=ref[:areas]
    Isles[i][:baseMVA]=ref[:baseMVA]
    Isles[i][:per_unit]=ref[:per_unit]
    Isles[i][:off_angmin]=ref[:off_angmin]
    Isles[i][:off_angmax]=ref[:off_angmax]


    #--
    Isles[i][:bus]=filter((x,y)->x in Isles[i][:map],ref[:bus])

    #--
    Isles[i][:gen]=filter((x,y)->y["gen_bus"] in Isles[i][:map],ref[:gen])

    #--
    Isles[i][:branch]=filter((x,y)->(y["f_bus"] in Isles[i][:map] || y["t_bus"] in Isles[i][:map]) && y["br_status"]==1   ,ref[:branch])

    #--
    Isles[i][:dcline]=filter((x,y)->(y["f_bus"] in Isles[i][:map] || y["t_bus"] in Isles[i][:map]) && y["br_status"]==1   ,ref[:dcline])    #-- untested



    #--
    Isles[i][:arcs_from] = [(jdx,branch["f_bus"],branch["t_bus"]) for (jdx,branch) in Isles[i][:branch]]
    Isles[i][:arcs_to]   = [(jdx,branch["t_bus"],branch["f_bus"]) for (jdx,branch) in Isles[i][:branch]]
    Isles[i][:arcs] = [Isles[i][:arcs_from]; Isles[i][:arcs_to]]

    #--
    Isles[i][:arcs_from_dc] = [(jdx,dcline["f_bus"],dcline["t_bus"]) for (jdx,dcline) in   Isles[i][:dcline]]
    Isles[i][:arcs_to_dc]   = [(jdx,dcline["t_bus"],dcline["f_bus"]) for (jdx,dcline) in   Isles[i][:dcline]]
    Isles[i][:arcs_dc]      = [  Isles[i][:arcs_from_dc];   Isles[i][:arcs_to_dc]]



    #--
    arcs_dc_param = Isles[i][:arcs_dc_param] = Dict()
    for (ldx,idx,jdx) in Isles[i][:arcs_from_dc]
      arcs_dc_param[(ldx,idx,jdx)] = Dict{String,Any}(
      "pmin" => Isles[i][:dcline][ldx]["pminf"],
      "pmax" => Isles[i][:dcline][ldx]["pmaxf"],
      "pref" => Isles[i][:dcline][ldx]["pf"],
      "qmin" => Isles[i][:dcline][ldx]["qminf"],
      "qmax" => Isles[i][:dcline][ldx]["qmaxf"],
      "qref" => Isles[i][:dcline][ldx]["qf"]
      )
      arcs_dc_param[(ldx,jdx,idx)] = Dict{String,Any}(
      "pmin" => Isles[i][:dcline][ldx]["pmint"],
      "pmax" => Isles[i][:dcline][ldx]["pmaxt"],
      "pref" => Isles[i][:dcline][ldx]["pt"],
      "qmin" => Isles[i][:dcline][ldx]["qmint"],
      "qmax" => Isles[i][:dcline][ldx]["qmaxt"],
      "qref" => Isles[i][:dcline][ldx]["qt"]
      )
    end


    #--
    bus_gens = Dict([(jdx, []) for (jdx,bus) in Isles[i][:bus]])
    for (jdx,gen) in Isles[i][:gen]
      push!(bus_gens[gen["gen_bus"]], jdx)
    end
    Isles[i][:bus_gens] = bus_gens

    #--
    bus_arcs = Dict([(idx, []) for (idx,bus) in Isles[i][:bus]])
    for (ldx,idx,jdx) in Isles[i][:arcs]
      push!(bus_arcs[idx], (ldx,idx,jdx))
    end
    Isles[i][:bus_arcs] = bus_arcs


    #--
    bus_arcs_dc = Dict([(idx, []) for (idx,bus) in Isles[i][:bus]])
    for (ldx,idx,jdx) in Isles[i][:arcs_dc]
      push!(bus_arcs_dc[idx], (ldx,idx,jdx))
    end
    Isles[i][:bus_arcs_dc] = bus_arcs_dc


    #--
    Isles[i][:buspairs] = buspair_parameters(Isles[i][:arcs_from], Isles[i][:branch], Isles[i][:bus])

    #--
    if length(keys(Isles[i][:gen]))>0

      #--
      ref_buses = Dict()
      big_gen = biggest_generator(Isles[i][:gen])
      gen_bus = big_gen["gen_bus"]
      ref_buses[gen_bus] = Isles[i][:bus][gen_bus]
      Isles[i][:ref_buses] = ref_buses


    else


      the_key=0


      for kdx in Isles[i][:map]
        the_key=kdx

        if the_key>0
          break
        end
      end
      Isles[i][:ref_buses]= filter((x,y)->y["bus_i"] ==the_key,Isles[i][:bus])
    end


  end


  #-----------#
  return Isles


end #-- function
