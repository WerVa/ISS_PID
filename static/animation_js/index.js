$(document).ready(function () {
	$('#retrieve').click(function () {
		let roller = document.querySelector('.lds-roller__wrapper');
		roller.classList.remove('hidden');
		$.ajax({
			url: "/Generate_plot",
			type: "GET",
			success: function (response) {
				$("#myimg").attr('src', '/static/' + response);
				roller.classList.add('hidden');
			},
			error: function (xhr) {
				//Do Something to handle error
			}
		});
	});

	$('#fuzzy').click(function () {
		let roller = document.querySelector('.lds-roller__wrapper');
		roller.classList.remove('hidden');
		$.ajax({
			url: "/Generate_plot_Fuzz",
			type: "GET",
			success: function (response) {
				$("#myimg").attr('src', '/static/' + response);
				roller.classList.add('hidden');
			},
			error: function (xhr) {
				//Do Something to handle error
			}
		});
	});


	$('#PidAndFuzzy').click(function () {
		let roller = document.querySelector('.lds-roller__wrapper');
		roller.classList.remove('hidden');
		$.ajax({
			url: "/Generate_porownanie",
			type: "GET",
			success: function (response) {
				$("#myimg").attr('src', '/static/' + response);
				roller.classList.add('hidden');
			},
			error: function (xhr) {
				//Do Something to handle error
			}
		});
	});
});
