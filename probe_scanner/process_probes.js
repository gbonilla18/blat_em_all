$(document).ready(function(){
    $('input:checkbox').attr('checked', true); 
	
	$('input:checkbox').each(function(){
		/*alert($(this).siblings('.multiple').length);*/
		
		if($(this).siblings('.multiple').length>0){
			$(this).attr('checked', false); 
			
		}
		
	});
	$('.download').click(function(){
		$('#probe_list').submit()
		
	});
});