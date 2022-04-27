#!/usr/bin/env bash

# --- colors --- #

#export NEWT_COLORS='
#root=gray,gray
#shadow=gray,gray
#window=,gray
#border=white,gray
#textbox=white,red
#button=black,white
#'

# --- functions --- #

m=0;
L=0;


main_menu() {

	whiptail --title "Main menu" --yesno --yes-button "Continue" --no-button "Exit" "Bienvenue dans l'assistant de résolution de \
	 l'équation de Shrodinger dans le cas du puit de potentiel. \n \nVous allez devoir choisir les paramètres initiaux" 15 40
	if [ $? -eq 0 ] ; then
		ask_params_1 
	else
		exit 
	fi
}

nan() {

	whiptail --title "error" --yesno --yes-button "Continue" --no-button "" "Veuillez entrer un nombre valide" 15 40
	if [[ $1 == 1 ]]; then
		ask_params_1
	else 
		ask_params_2
	fi
}

ask_params_1() {

	CHOICE=$(whiptail --title "1/2" --inputbox  "masse de la particule (eV)" --cancel-button "Back" 8 40 \
	3>&1 1>&2 2>&3)

	choice=${CHOICE,,}

	if [ $? -eq 0 ] ; then
		re='^[0-9]+$'
		if ! [[ $choice =~ $re ]] ; then
   			main_menu
		else
			$m = $choice
			ask_params_2
		fi
		 
	else
		main_menu
	fi

}

ask_params_2() {

	CHOICE=$(whiptail --title "2/2" --inputbox  "Distance entre les murs (nm)" --cancel-button "Back" 8 40 \
	3>&1 1>&2 2>&3)

	if [ $? -eq 0 ] ; then
		$L = ${CHOICE,,} 
	else
		ask_params_1 
	fi

}

print_graph() {
	i=1
	cat /dev/null > data_set.dat
	cat /dev/null > temp.p

	if [ "$mode" == "country" ] ; then

		while [ -n "${Country_Consumption[$i,$index]}" ] ; do
			echo -e "${Country_Consumption[$i,0]} \t ${Country_Consumption[$i,$index]}" >> data_set.dat
			((i++))
		done

	else

		while [ -n "${Continent_Consumption[$i,$index]}" ] ; do
			echo -e "${Continent_Consumption[$i,0]} \t ${Continent_Consumption[$i,$index]}" >> data_set.dat
			((i++))
		done

	fi

	echo -e "set terminal png size 1920,1080 \n set output '$choice.png' \n set title 'Consumption of "$CHOICE" by years (TWH)' \n plot \"data_set.dat\" w lp " >> temp.p

	mkdir -p results/"$choice" 

	gnuplot temp.p

	mv "$choice.png" results/"$choice"

	rm temp.p data_set.dat

	end
}

info_menu() {
	CHOICE=$(whiptail --title "Specific informations" --menu --notags "Want to know about : " --cancel-button "Back" 11 60 4 \
	"1" "The country who produce the more \"green\" electricity"  \
	"2" "The country who produce the more \"fossil\" electricity"  \
	"3" "The production of renewable energy"  \
	"4" "The comparison of renewable / fossil energy"  \
	3>&1 1>&2 2>&3)

	if [ $? -eq 1 ] ; then
		main_menu
	fi

	case $CHOICE in
		"1" ) green_max ;;
		"2" ) fossil_max ;;
		"3" ) renewable_power_gen_per_year ;;
		"4" ) global_generation ;;
	esac 
}

main_menu
